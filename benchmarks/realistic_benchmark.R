## realistic_benchmark.R — more realistic simulation + benchmark
## Exploratory: does RMTnet beat hard thresholding when the data isn't clean?

user_lib <- file.path(Sys.getenv("USERPROFILE"), "R", "library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

library(RMThreshold)
library(mclust)
library(dynamicTreeCut)

pkg_path <- file.path(dirname(getwd()), "test", "RMTnet")
for (f in list.files(file.path(pkg_path, "R"), full.names = TRUE)) source(f)

# ══════════════════════════════════════════════════════════════════════════════
# REALISTIC SIMULATOR
# ══════════════════════════════════════════════════════════════════════════════

simulate_realistic <- function(
    n_genes       = 500,
    n_samples     = 80,
    module_sizes  = c(60, 45, 35, 25, 15),    # unequal sizes
    overlap_frac  = 0.1,     # 10% of module genes belong to 2 modules
    hub_genes     = 5,       # genes that load on ALL modules
    loading_sd    = 0.25,    # within-module loading varies (mean=1, sd=0.25)
    bg_strength   = 0.15,    # background correlation (batch/housekeeping)
    signal_range  = c(0.4, 0.8),  # modules vary in signal strength
    seed          = 42
) {
  set.seed(seed)
  n_modules  <- length(module_sizes)

  # True labels: 0 = unassigned
  true_labels <- rep(0L, n_genes)

  # Assign genes to modules (non-overlapping core)
  idx <- 1
  module_gene_lists <- list()
  for (m in seq_len(n_modules)) {
    end <- min(idx + module_sizes[m] - 1, n_genes)
    module_gene_lists[[m]] <- idx:end
    true_labels[idx:end]   <- m
    idx <- end + 1
  }

  # Add overlapping genes: some genes from module m also load on module m+1
  for (m in seq_len(n_modules - 1)) {
    n_overlap <- max(1, round(module_sizes[m] * overlap_frac))
    overlap_genes <- sample(module_gene_lists[[m]], n_overlap)
    module_gene_lists[[m + 1]] <- c(module_gene_lists[[m + 1]], overlap_genes)
    # These genes keep their primary label but also participate in next module
  }

  # Hub genes: random genes that load weakly on every module
  hub_idx <- sample(which(true_labels == 0),
                    min(hub_genes, sum(true_labels == 0)))

  # Generate latent factors (one per module)
  latent <- matrix(rnorm(n_modules * n_samples), nrow = n_modules)

  # Background factor (batch effect / housekeeping)
  bg_factor <- rnorm(n_samples)

  # Generate expression
  mat <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)

  # Per-module signal strengths (vary between modules)
  mod_strengths <- seq(signal_range[1], signal_range[2], length.out = n_modules)

  for (m in seq_len(n_modules)) {
    genes_in_mod <- module_gene_lists[[m]]
    s <- mod_strengths[m]

    for (g in genes_in_mod) {
      # Variable loading: each gene responds differently to the latent factor
      loading <- abs(rnorm(1, mean = 1, sd = loading_sd))
      mat[g, ] <- sqrt(1 - s) * mat[g, ] + sqrt(s) * loading * latent[m, ]
    }
  }

  # Hub genes: weak loading on every module
  for (g in hub_idx) {
    hub_load <- 0.3  # weaker than module-specific genes
    combined <- colMeans(latent)  # average of all latent factors
    mat[g, ] <- sqrt(1 - hub_load) * mat[g, ] + sqrt(hub_load) * combined
  }

  # Background correlation: every gene gets a small batch effect
  for (g in seq_len(n_genes)) {
    mat[g, ] <- mat[g, ] + bg_strength * bg_factor
  }

  # Standardise rows
  mat <- t(scale(t(mat)))

  rownames(mat) <- paste0("gene", seq_len(n_genes))
  colnames(mat) <- paste0("sample", seq_len(n_samples))

  list(mat = mat, true_labels = true_labels, n_modules = n_modules,
       module_sizes = module_sizes, hub_idx = hub_idx)
}

# ══════════════════════════════════════════════════════════════════════════════
# BENCHMARK HELPERS
# ══════════════════════════════════════════════════════════════════════════════

cluster_corr <- function(corr, min_module_size = 15) {
  n     <- nrow(corr)
  d_mat <- matrix(sqrt(pmax(0, 2 * (1 - corr))), n, n)
  hc    <- hclust(as.dist(d_mat), method = "average")
  cutreeDynamic(hc, distM = d_mat, minClusterSize = min_module_size, verbose = 0)
}

auto_rmth_threshold <- function(res) {
  valid <- which(res$p.ks > 0.05)
  if (length(valid) > 0) return(res$tested.thresholds[min(valid)])
  res$tested.thresholds[which.min(res$sse.exp)]
}

run_trial <- function(sim, min_mod = 15) {
  mat   <- sim$mat
  truth <- sim$true_labels
  result <- list()

  # 1. RMTnet (spectral filtering)
  tryCatch({
    net  <- rmt_network(mat, min_module_size = min_mod, verbose = FALSE)
    result$rmtnet <- adjustedRandIndex(truth, net$modules)
  }, error = function(e) { result$rmtnet <<- NA_real_ })

  # 2. RMThreshold (hard threshold via NNSD)
  tryCatch({
    corr_raw <- cor(t(mat))
    res      <- rm.get.threshold(corr_raw, interactive = FALSE,
                                 plot.comp = FALSE, save.fit = FALSE,
                                 plot.spacing = FALSE)
    thresh   <- auto_rmth_threshold(res)
    corr_th  <- corr_raw
    corr_th[abs(corr_th) < thresh] <- 0
    diag(corr_th) <- 1
    result$rmthreshold <- adjustedRandIndex(truth, cluster_corr(corr_th, min_mod))
  }, error = function(e) { result$rmthreshold <<- NA_real_ })

  # 3. Baseline (raw correlation, no filtering)
  tryCatch({
    corr_raw <- cor(t(mat))
    result$baseline <- adjustedRandIndex(truth, cluster_corr(corr_raw, min_mod))
  }, error = function(e) { result$baseline <<- NA_real_ })

  result
}

# ══════════════════════════════════════════════════════════════════════════════
# RUN BENCHMARK
# ══════════════════════════════════════════════════════════════════════════════

out_file <- "C:/Users/felix/downloads/test/realistic_results.txt"
sink(out_file, split = TRUE)

cat("Realistic Simulation Benchmark\n")
cat("Features: overlapping modules, variable loadings, unequal module\n")
cat("          sizes, hub genes, weak signal regimes\n\n")

scenarios <- list(
  list(label = "Messy (all features)",
       n_genes = 500, n_samples = 80,
       module_sizes = c(60, 45, 35, 25, 15),
       overlap_frac = 0.10, hub_genes = 5,
       loading_sd = 0.25, bg_strength = 0.0, signal_range = c(0.4, 0.8)),

  list(label = "Heavy overlap (20%)",
       n_genes = 500, n_samples = 80,
       module_sizes = c(60, 45, 35, 25, 15),
       overlap_frac = 0.20, hub_genes = 8,
       loading_sd = 0.30, bg_strength = 0.0, signal_range = c(0.4, 0.7)),

  list(label = "Many hub genes (15)",
       n_genes = 500, n_samples = 80,
       module_sizes = c(60, 45, 35, 25, 15),
       overlap_frac = 0.10, hub_genes = 15,
       loading_sd = 0.25, bg_strength = 0.0, signal_range = c(0.4, 0.8)),

  list(label = "Many small modules",
       n_genes = 500, n_samples = 100,
       module_sizes = c(30, 25, 20, 20, 15, 15, 15, 10),
       overlap_frac = 0.10, hub_genes = 10,
       loading_sd = 0.20, bg_strength = 0.0, signal_range = c(0.5, 0.8)),

  list(label = "Weak signal + noise",
       n_genes = 500, n_samples = 60,
       module_sizes = c(50, 40, 30, 20),
       overlap_frac = 0.15, hub_genes = 5,
       loading_sd = 0.30, bg_strength = 0.0, signal_range = c(0.25, 0.5))
)

N_REPS <- 10
all_results <- list()

for (sc in scenarios) {
  cat(sprintf("Scenario: %s\n", sc$label))

  reps <- lapply(seq_len(N_REPS), function(i) {
    sim <- simulate_realistic(
      n_genes       = sc$n_genes,
      n_samples     = sc$n_samples,
      module_sizes  = sc$module_sizes,
      overlap_frac  = sc$overlap_frac,
      hub_genes     = sc$hub_genes,
      loading_sd    = sc$loading_sd,
      bg_strength   = sc$bg_strength,
      signal_range  = sc$signal_range,
      seed          = i * 100
    )
    run_trial(sim, min_mod = 10)
  })

  rmtnet_ari   <- sapply(reps, `[[`, "rmtnet")
  rmthresh_ari <- sapply(reps, `[[`, "rmthreshold")
  baseline_ari <- sapply(reps, `[[`, "baseline")

  cat(sprintf("  %-18s  ARI = %.3f +/- %.3f\n",
              "RMTnet",      mean(rmtnet_ari,   na.rm=TRUE), sd(rmtnet_ari,   na.rm=TRUE)))
  cat(sprintf("  %-18s  ARI = %.3f +/- %.3f\n",
              "RMThreshold", mean(rmthresh_ari, na.rm=TRUE), sd(rmthresh_ari, na.rm=TRUE)))
  cat(sprintf("  %-18s  ARI = %.3f +/- %.3f\n\n",
              "Baseline",    mean(baseline_ari, na.rm=TRUE), sd(baseline_ari, na.rm=TRUE)))

  all_results[[sc$label]] <- list(
    rmtnet = rmtnet_ari, rmthreshold = rmthresh_ari, baseline = baseline_ari
  )
}

# Summary
all_rmt  <- unlist(lapply(all_results, `[[`, "rmtnet"))
all_rmth <- unlist(lapply(all_results, `[[`, "rmthreshold"))
all_base <- unlist(lapply(all_results, `[[`, "baseline"))

cat("══════════════════════════════════════════════════════════════════\n")
cat("OVERALL\n")
cat(sprintf("  RMTnet       mean ARI = %.3f\n", mean(all_rmt,  na.rm=TRUE)))
cat(sprintf("  RMThreshold  mean ARI = %.3f\n", mean(all_rmth, na.rm=TRUE)))
cat(sprintf("  Baseline     mean ARI = %.3f\n", mean(all_base, na.rm=TRUE)))
cat("══════════════════════════════════════════════════════════════════\n")
cat("\nARI: 1.0 = perfect, 0 = random, <0 = worse than random\n")

# Comparison with original clean simulation
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("COMPARISON: Clean vs Realistic simulation\n")
cat("──────────────────────────────────────────────────────────────────\n")
cat("                    Clean sim    Realistic sim\n")
cat(sprintf("  RMTnet              0.302         %.3f\n", mean(all_rmt,  na.rm=TRUE)))
cat(sprintf("  RMThreshold         0.904         %.3f\n", mean(all_rmth, na.rm=TRUE)))
cat(sprintf("  Baseline            0.206         %.3f\n", mean(all_base, na.rm=TRUE)))
cat("══════════════════════════════════════════════════════════════════\n")

sink()
cat("\nResults saved to", out_file, "\n")

# Cleanup RMThreshold PNGs
invisible(file.remove(list.files(pattern = "^RMT.*\\.png$")))
invisible(file.remove(list.files(pattern = "^Rplots\\.pdf$")))
