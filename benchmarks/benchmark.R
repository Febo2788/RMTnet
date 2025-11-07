## benchmark.R — compare RMTnet vs RMThreshold vs Baseline
## Metric: Adjusted Rand Index (ARI) — 1.0 = perfect module recovery, 0 = random

user_lib <- file.path(Sys.getenv("USERPROFILE"), "R", "library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

for (p in c("RMThreshold", "mclust", "dynamicTreeCut")) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, lib = user_lib, repos = "https://cloud.r-project.org", quiet = TRUE)
}

pkg_path <- file.path(dirname(getwd()), "test", "RMTnet")
for (f in list.files(file.path(pkg_path, "R"), full.names = TRUE)) source(f)

library(RMThreshold)
library(mclust)
library(dynamicTreeCut)

# ── helpers ──────────────────────────────────────────────────────────────────

ari <- function(true_labels, pred_labels) adjustedRandIndex(true_labels, pred_labels)

true_labels <- function(n_genes, n_modules, module_size) {
  labels <- rep(0L, n_genes)
  for (m in seq_len(n_modules)) {
    idx <- ((m - 1) * module_size + 1):(m * module_size)
    labels[idx] <- m
  }
  labels
}

cluster_corr <- function(corr, min_module_size = 20) {
  n     <- nrow(corr)
  d_mat <- matrix(sqrt(pmax(0, 2 * (1 - corr))), n, n)
  hc    <- hclust(as.dist(d_mat), method = "average")
  cutreeDynamic(hc, distM = d_mat, minClusterSize = min_module_size, verbose = 0)
}

# Auto-select RMThreshold threshold:
# find the first threshold where p_ks > 0.05 (onset of GOE-like spacing)
auto_rmth_threshold <- function(res) {
  valid <- which(res$p.ks > 0.05)
  if (length(valid) > 0)
    return(res$tested.thresholds[min(valid)])
  # fallback: minimum SSE from exponential fit
  res$tested.thresholds[which.min(res$sse.exp)]
}

# ── single trial ─────────────────────────────────────────────────────────────

run_trial <- function(n_genes, n_samples, n_modules, module_size,
                      signal_strength, seed, min_mod = 20) {
  mat   <- simulate_expression(n_genes, n_samples, n_modules,
                               module_size, signal_strength, seed)
  truth <- true_labels(n_genes, n_modules, module_size)
  result <- list()

  # ── 1. RMTnet (spectral filtering) ──────────────────────────────────────
  tryCatch({
    net  <- rmt_network(mat, min_module_size = min_mod, verbose = FALSE)
    result$rmtnet <- ari(truth, net$modules)
  }, error = function(e) { result$rmtnet <<- NA_real_ })

  # ── 2. RMThreshold (hard threshold via NNSD / GOE transition) ───────────
  tryCatch({
    corr_raw <- cor(t(mat))
    # suppress PNG output and interactivity
    res    <- rm.get.threshold(corr_raw, interactive = FALSE,
                               plot.comp = FALSE, save.fit = FALSE,
                               plot.spacing = FALSE)
    thresh <- auto_rmth_threshold(res)
    corr_th <- corr_raw
    corr_th[abs(corr_th) < thresh] <- 0
    diag(corr_th) <- 1
    result$rmthreshold <- ari(truth, cluster_corr(corr_th, min_mod))
  }, error = function(e) { result$rmthreshold <<- NA_real_ })

  # ── 3. Baseline (raw correlation, no noise removal) ──────────────────────
  tryCatch({
    corr_raw <- cor(t(mat))
    result$baseline <- ari(truth, cluster_corr(corr_raw, min_mod))
  }, error = function(e) { result$baseline <<- NA_real_ })

  result
}

# ── experiment grid ───────────────────────────────────────────────────────────

cat("Running benchmarks (10 reps x 5 scenarios)...\n\n")

scenarios <- list(
  list(label = "Easy   (high signal)",   n_genes=300, n_samples=80,  n_modules=3, module_size=50, signal=0.8),
  list(label = "Medium (typical)   ",    n_genes=400, n_samples=80,  n_modules=4, module_size=40, signal=0.6),
  list(label = "Hard   (low signal)",    n_genes=500, n_samples=60,  n_modules=4, module_size=40, signal=0.4),
  list(label = "Small N (few samples)",  n_genes=400, n_samples=40,  n_modules=3, module_size=40, signal=0.6),
  list(label = "Many modules       ",    n_genes=500, n_samples=100, n_modules=6, module_size=40, signal=0.7)
)

N_REPS <- 10

all_results <- list()

for (sc in scenarios) {
  cat(sprintf("Scenario: %s\n", sc$label))
  reps <- lapply(seq_len(N_REPS), function(i) {
    run_trial(sc$n_genes, sc$n_samples, sc$n_modules, sc$module_size,
              sc$signal, seed = i * 100, min_mod = 20)
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

  all_results[[sc$label]] <- list(rmtnet=rmtnet_ari, rmthreshold=rmthresh_ari, baseline=baseline_ari)
}

# ── overall summary ───────────────────────────────────────────────────────────

all_rmt  <- unlist(lapply(all_results, `[[`, "rmtnet"))
all_rmth <- unlist(lapply(all_results, `[[`, "rmthreshold"))
all_base <- unlist(lapply(all_results, `[[`, "baseline"))

cat("══════════════════════════════════════════\n")
cat("OVERALL (all scenarios, all reps)\n")
cat(sprintf("  RMTnet       mean ARI = %.3f\n", mean(all_rmt,  na.rm=TRUE)))
cat(sprintf("  RMThreshold  mean ARI = %.3f\n", mean(all_rmth, na.rm=TRUE)))
cat(sprintf("  Baseline     mean ARI = %.3f\n", mean(all_base, na.rm=TRUE)))
cat("══════════════════════════════════════════\n")
cat("\nARI: 1.0 = perfect, 0 = no better than random, <0 = worse than random\n")

# cleanup PNGs dropped by RMThreshold
invisible(file.remove(list.files(pattern = "^RMT.*\\.png$")))
