## fuzzy_benchmark.R — continuous module membership + correlation recovery metric

user_lib <- file.path(Sys.getenv("USERPROFILE"), "R", "library")
.libPaths(c(user_lib, .libPaths()))
library(RMThreshold)
library(dynamicTreeCut)

pkg_path <- file.path(dirname(getwd()), "test", "RMTnet")
for (f in list.files(file.path(pkg_path, "R"), full.names = TRUE)) source(f)

# ══════════════════════════════════════════════════════════════════════════════
# FUZZY SIMULATOR
# ══════════════════════════════════════════════════════════════════════════════
#
# Instead of: gene belongs to module A (1.0) or not (0.0)
# Now:        gene loads 0.6 on module A, 0.3 on module B, 0.1 on module C
#
# This creates continuous correlation gradients — no clean block boundaries.

simulate_fuzzy <- function(n_genes = 500, n_samples = 80, n_factors = 5,
                           sparsity = 0.6, seed = 42) {
  set.seed(seed)

  # Latent factors
  latent <- matrix(rnorm(n_factors * n_samples), nrow = n_factors)

  # Loading matrix: each gene loads on multiple factors
  # Dirichlet-like: sample from exponential, normalize, then sparsify
  W <- matrix(rexp(n_genes * n_factors, rate = 1), nrow = n_genes)

  # Sparsify: randomly zero out some loadings so genes don't load on everything
  mask <- matrix(rbinom(n_genes * n_factors, 1, prob = 1 - sparsity),
                 nrow = n_genes)
  # Ensure every gene loads on at least one factor
  for (g in seq_len(n_genes)) {
    if (sum(mask[g, ]) == 0) mask[g, sample(n_factors, 1)] <- 1
  }
  W <- W * mask

  # Normalize rows to sum to 1 (relative contribution of each factor)
  W <- W / rowSums(W)

  # Signal strength varies per gene (some genes are noisier)
  signal_strength <- runif(n_genes, 0.3, 0.8)

  # Generate expression: gene = signal * (weighted sum of factors) + noise
  signal_part <- W %*% latent  # n_genes x n_samples
  noise_part  <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)

  mat <- sweep(signal_part, 1, sqrt(signal_strength), "*") +
         sweep(noise_part,  1, sqrt(1 - signal_strength), "*")
  mat <- t(scale(t(mat)))

  # TRUE correlation: what the correlation matrix would look like with
  # infinite samples (no noise). This is the signal we're trying to recover.
  # C_true = diag(s) * W * W^T * diag(s) + diag(1-s)
  # But normalized to correlation scale:
  S      <- diag(signal_strength)
  C_true <- S %*% tcrossprod(W) %*% S
  diag(C_true) <- 1
  # Ensure it's a proper correlation matrix
  D_inv <- diag(1 / sqrt(diag(C_true)))
  C_true <- D_inv %*% C_true %*% D_inv

  rownames(mat) <- paste0("gene", seq_len(n_genes))
  colnames(mat) <- paste0("sample", seq_len(n_samples))

  list(mat = mat, W = W, C_true = C_true, signal_strength = signal_strength)
}

# ══════════════════════════════════════════════════════════════════════════════
# METRIC: correlation matrix recovery
# ══════════════════════════════════════════════════════════════════════════════
#
# Frobenius norm of the difference, normalized:
# recovery_error = ||C_method - C_true||_F / ||C_true||_F
# Lower = better. 0 = perfect recovery.

recovery_error <- function(C_method, C_true) {
  norm(C_method - C_true, type = "F") / norm(C_true, type = "F")
}

# Also measure: correlation between upper-triangle entries
# (do the relative strengths match?)
recovery_cor <- function(C_method, C_true) {
  idx <- upper.tri(C_true)
  cor(C_method[idx], C_true[idx])
}

# ══════════════════════════════════════════════════════════════════════════════
# HELPER: auto threshold for RMThreshold
# ══════════════════════════════════════════════════════════════════════════════

auto_rmth_threshold <- function(res) {
  valid <- which(res$p.ks > 0.05)
  if (length(valid) > 0) return(res$tested.thresholds[min(valid)])
  res$tested.thresholds[which.min(res$sse.exp)]
}

# ══════════════════════════════════════════════════════════════════════════════
# SINGLE TRIAL
# ══════════════════════════════════════════════════════════════════════════════

run_trial <- function(sim) {
  mat    <- sim$mat
  C_true <- sim$C_true
  result <- list()

  # Raw (no filtering)
  C_raw <- cor(t(mat))
  result$raw_error <- recovery_error(C_raw, C_true)
  result$raw_cor   <- recovery_cor(C_raw, C_true)

  # RMTnet (spectral filtering)
  tryCatch({
    mp   <- mp_threshold(mat)
    filt <- spectral_filter(mp)
    result$rmtnet_error <- recovery_error(filt$cleaned_corr, C_true)
    result$rmtnet_cor   <- recovery_cor(filt$cleaned_corr, C_true)
  }, error = function(e) {
    result$rmtnet_error <<- NA; result$rmtnet_cor <<- NA
  })

  # RMThreshold (hard threshold)
  tryCatch({
    res    <- rm.get.threshold(C_raw, interactive = FALSE,
                               plot.comp = FALSE, save.fit = FALSE,
                               plot.spacing = FALSE)
    thresh <- auto_rmth_threshold(res)
    C_hard <- C_raw
    C_hard[abs(C_hard) < thresh] <- 0
    diag(C_hard) <- 1
    result$rmth_error <- recovery_error(C_hard, C_true)
    result$rmth_cor   <- recovery_cor(C_hard, C_true)
  }, error = function(e) {
    result$rmth_error <<- NA; result$rmth_cor <<- NA
  })

  result
}

# ══════════════════════════════════════════════════════════════════════════════
# RUN BENCHMARK
# ══════════════════════════════════════════════════════════════════════════════

out_file <- "C:/Users/felix/downloads/test/fuzzy_results.txt"
sink(out_file, split = TRUE)

cat("Fuzzy Module Benchmark — Correlation Matrix Recovery\n")
cat("Genes load continuously on multiple factors (no clean blocks).\n")
cat("Metric: how close is the method's output to the true signal matrix?\n\n")

scenarios <- list(
  list(label = "Standard fuzzy (5 factors)",
       n_genes = 500, n_samples = 80, n_factors = 5, sparsity = 0.6),
  list(label = "Dense loadings (low sparsity)",
       n_genes = 500, n_samples = 80, n_factors = 5, sparsity = 0.3),
  list(label = "Many factors (10)",
       n_genes = 500, n_samples = 100, n_factors = 10, sparsity = 0.7),
  list(label = "Few samples (Q=0.08)",
       n_genes = 500, n_samples = 40, n_factors = 5, sparsity = 0.6),
  list(label = "Large + sparse (1000 genes)",
       n_genes = 1000, n_samples = 100, n_factors = 8, sparsity = 0.7)
)

N_REPS <- 10

for (sc in scenarios) {
  cat(sprintf("Scenario: %s\n", sc$label))

  reps <- lapply(seq_len(N_REPS), function(i) {
    sim <- simulate_fuzzy(sc$n_genes, sc$n_samples, sc$n_factors,
                          sc$sparsity, seed = i * 100)
    run_trial(sim)
  })

  # Recovery error (lower = better)
  raw_err  <- sapply(reps, `[[`, "raw_error")
  rmt_err  <- sapply(reps, `[[`, "rmtnet_error")
  rmth_err <- sapply(reps, `[[`, "rmth_error")

  # Recovery correlation (higher = better)
  raw_cor  <- sapply(reps, `[[`, "raw_cor")
  rmt_cor  <- sapply(reps, `[[`, "rmtnet_cor")
  rmth_cor <- sapply(reps, `[[`, "rmth_cor")

  cat("  Recovery error (lower = better):\n")
  cat(sprintf("    %-18s  %.4f +/- %.4f\n", "RMTnet",      mean(rmt_err,  na.rm=TRUE), sd(rmt_err,  na.rm=TRUE)))
  cat(sprintf("    %-18s  %.4f +/- %.4f\n", "RMThreshold", mean(rmth_err, na.rm=TRUE), sd(rmth_err, na.rm=TRUE)))
  cat(sprintf("    %-18s  %.4f +/- %.4f\n", "Raw (no filt)", mean(raw_err, na.rm=TRUE), sd(raw_err, na.rm=TRUE)))

  cat("  Recovery correlation (higher = better):\n")
  cat(sprintf("    %-18s  %.4f +/- %.4f\n", "RMTnet",       mean(rmt_cor,  na.rm=TRUE), sd(rmt_cor,  na.rm=TRUE)))
  cat(sprintf("    %-18s  %.4f +/- %.4f\n", "RMThreshold",  mean(rmth_cor, na.rm=TRUE), sd(rmth_cor, na.rm=TRUE)))
  cat(sprintf("    %-18s  %.4f +/- %.4f\n\n", "Raw (no filt)", mean(raw_cor, na.rm=TRUE), sd(raw_cor, na.rm=TRUE)))
}

cat("══════════════════════════════════════════════════════════════════\n")
cat("INTERPRETATION\n")
cat("──────────────────────────────────────────────────────────────────\n")
cat("Recovery error:  how far is the method's matrix from the truth\n")
cat("                 (Frobenius norm, normalized). Lower = better.\n")
cat("Recovery corr:   do the pairwise correlations in the method's\n")
cat("                 output track the true correlations? Higher = better.\n")
cat("══════════════════════════════════════════════════════════════════\n")

sink()
cat("\nResults saved to", out_file, "\n")

invisible(file.remove(list.files(pattern = "^RMT.*\\.png$")))
invisible(file.remove(list.files(pattern = "^Rplots\\.pdf$")))
