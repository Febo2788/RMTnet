## batch_explainer.R — visual explanation of why batch effects break spectral filtering

user_lib <- file.path(Sys.getenv("USERPROFILE"), "R", "library")
.libPaths(c(user_lib, .libPaths()))
library(dynamicTreeCut)
library(mclust)

pkg_path <- file.path(dirname(getwd()), "test", "RMTnet")
for (f in list.files(file.path(pkg_path, "R"), full.names = TRUE)) source(f)

out_dir <- "C:/Users/felix/downloads/test"

# ═══════════════════════════════════════════════════════════════════════════
# Generate three versions of the same data
# ═══════════════════════════════════════════════════════════════════════════

set.seed(42)
n_genes   <- 300
n_samples <- 80
n_modules <- 3
mod_sizes <- c(60, 60, 60)

# Latent factors (one per module)
latent <- matrix(rnorm(n_modules * n_samples), nrow = n_modules)

# Batch factor — the villain
batch_factor <- rnorm(n_samples)

# Base noise for every gene (same across all 3 versions)
noise <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)

build_mat <- function(noise, latent, batch_strength) {
  mat <- noise
  signal <- 0.7
  # Module 1: genes 1-60
  for (g in 1:60)
    mat[g, ] <- sqrt(1-signal)*mat[g, ] + sqrt(signal)*latent[1, ]
  # Module 2: genes 61-120
  for (g in 61:120)
    mat[g, ] <- sqrt(1-signal)*mat[g, ] + sqrt(signal)*latent[2, ]
  # Module 3: genes 121-180
  for (g in 121:180)
    mat[g, ] <- sqrt(1-signal)*mat[g, ] + sqrt(signal)*latent[3, ]
  # Add batch to ALL genes
  for (g in 1:n_genes)
    mat[g, ] <- mat[g, ] + batch_strength * batch_factor
  t(scale(t(mat)))
}

mat_clean  <- build_mat(noise, latent, batch_strength = 0.00)
mat_mild   <- build_mat(noise, latent, batch_strength = 0.15)
mat_strong <- build_mat(noise, latent, batch_strength = 0.35)

true_labels <- c(rep(1,60), rep(2,60), rep(3,60), rep(0, 120))

# ═══════════════════════════════════════════════════════════════════════════
# Run RMTnet on each
# ═══════════════════════════════════════════════════════════════════════════

mp_clean  <- mp_threshold(mat_clean)
mp_mild   <- mp_threshold(mat_mild)
mp_strong <- mp_threshold(mat_strong)

net_clean  <- rmt_network(mat_clean,  min_module_size = 15, verbose = FALSE)
net_mild   <- rmt_network(mat_mild,   min_module_size = 15, verbose = FALSE)
net_strong <- rmt_network(mat_strong, min_module_size = 15, verbose = FALSE)

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 1: Correlation heatmaps — what the raw data looks like
# ═══════════════════════════════════════════════════════════════════════════

png(file.path(out_dir, "fig1_correlation_heatmaps.png"),
    width = 1400, height = 500, res = 130)
par(mfrow = c(1, 3), mar = c(3, 3, 4, 2))

heatcol <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

cor_clean  <- cor(t(mat_clean))
cor_mild   <- cor(t(mat_mild))
cor_strong <- cor(t(mat_strong))

image(cor_clean[, nrow(cor_clean):1], col = heatcol, zlim = c(-1,1),
      axes = FALSE, main = "No batch effect")
mtext("Modules are visible as red blocks on the diagonal.\nBackground is white (zero correlation).",
      side = 1, line = 1.5, cex = 0.65, col = "grey30")

image(cor_mild[, nrow(cor_mild):1], col = heatcol, zlim = c(-1,1),
      axes = FALSE, main = "Mild batch (0.15)")
mtext("Blocks still visible, but the background has turned\nslightly pink — everything is a little correlated now.",
      side = 1, line = 1.5, cex = 0.65, col = "grey30")

image(cor_strong[, nrow(cor_strong):1], col = heatcol, zlim = c(-1,1),
      axes = FALSE, main = "Strong batch (0.35)")
mtext("The whole matrix is pink. The batch correlates every gene\nwith every other gene. Module blocks are harder to see.",
      side = 1, line = 1.5, cex = 0.65, col = "grey30")
dev.off()

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 2: Eigenvalue spectra — what the MP law sees
# ═══════════════════════════════════════════════════════════════════════════

png(file.path(out_dir, "fig2_eigenvalue_spectra.png"),
    width = 1400, height = 500, res = 130)
par(mfrow = c(1, 3), mar = c(4, 4, 4, 2))

plot_ev <- function(mp_res, title, subtitle) {
  ev   <- sort(mp_res$eigenvalues, decreasing = TRUE)
  lp   <- mp_res$lambda_plus
  cols <- ifelse(ev > lp, "#D62728", "#AAAAAA")
  plot(seq_along(ev), ev, pch = 16, cex = 0.7, col = cols,
       xlab = "PC rank", ylab = "Eigenvalue",
       main = title, log = "y", ylim = c(0.01, max(ev)*1.5))
  abline(h = lp, col = "red", lty = 2, lwd = 2)
  text(length(ev)*0.7, lp*1.5,
       paste0("lambda+ = ", round(lp, 3)), col = "red", cex = 0.8)
  n_sig <- sum(ev > lp)
  mtext(paste0(subtitle, "\n", n_sig, " signal PCs (red dots above the line)"),
        side = 1, line = 2.5, cex = 0.6, col = "grey30")
}

plot_ev(mp_clean,  "No batch effect",
        "3 big eigenvalues = 3 modules. Everything else is below the line (noise).")
plot_ev(mp_mild,   "Mild batch (0.15)",
        "Still 3 big ones, but the largest grew — it absorbed some batch signal.")
plot_ev(mp_strong, "Strong batch (0.35)",
        "One MASSIVE eigenvalue dominates. That's the batch, not biology.\nThe MP law can't tell the difference — it just sees 'statistically real'.")
dev.off()

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 3: What spectral filtering does to the correlation matrix
# ═══════════════════════════════════════════════════════════════════════════

filt_clean  <- spectral_filter(mp_clean)
filt_strong <- spectral_filter(mp_strong)

png(file.path(out_dir, "fig3_spectral_filter_result.png"),
    width = 1400, height = 500, res = 130)
par(mfrow = c(1, 3), mar = c(3, 3, 4, 2))

image(filt_clean$cleaned_corr[, n_genes:1], col = heatcol, zlim = c(-1,1),
      axes = FALSE, main = "Cleaned — no batch")
mtext("Spectral filtering removed the noise.\nModule blocks are crisp. This is the ideal case.",
      side = 1, line = 1.5, cex = 0.65, col = "grey30")

image(filt_strong$cleaned_corr[, n_genes:1], col = heatcol, zlim = c(-1,1),
      axes = FALSE, main = "Cleaned — strong batch")
mtext("Spectral filtering kept the batch because it looked real.\nThe whole matrix is still pink. Modules are buried.",
      side = 1, line = 1.5, cex = 0.65, col = "grey30")

# Show what it SHOULD look like: remove batch first, then filter
mat_corrected <- mat_strong
# Simple batch correction: regress out the batch factor from each gene
for (g in 1:n_genes) {
  fit <- lm(mat_corrected[g, ] ~ batch_factor)
  mat_corrected[g, ] <- residuals(fit)
}
mat_corrected <- t(scale(t(mat_corrected)))
mp_corrected   <- mp_threshold(mat_corrected)
filt_corrected <- spectral_filter(mp_corrected)

image(filt_corrected$cleaned_corr[, n_genes:1], col = heatcol, zlim = c(-1,1),
      axes = FALSE, main = "Remove batch FIRST, then filter")
mtext("If you regress out the batch before running RMTnet,\nthe modules come back. This is the correct workflow.",
      side = 1, line = 1.5, cex = 0.65, col = "grey30")
dev.off()

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 4: Module recovery comparison
# ═══════════════════════════════════════════════════════════════════════════

png(file.path(out_dir, "fig4_module_recovery.png"),
    width = 1000, height = 500, res = 130)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# RMTnet ARI across batch strengths
batch_vals <- seq(0, 0.5, by = 0.05)
ari_rmtnet <- numeric(length(batch_vals))
ari_rmth   <- numeric(length(batch_vals))

for (i in seq_along(batch_vals)) {
  m <- build_mat(noise, latent, batch_strength = batch_vals[i])

  # RMTnet
  tryCatch({
    n <- rmt_network(m, min_module_size = 15, verbose = FALSE)
    ari_rmtnet[i] <- adjustedRandIndex(true_labels, n$modules)
  }, error = function(e) ari_rmtnet[i] <<- NA)

  # Hard threshold at r=0.3 (simple stand-in for RMThreshold)
  tryCatch({
    cr <- cor(t(m))
    cr[abs(cr) < 0.3] <- 0
    diag(cr) <- 1
    ari_rmth[i] <- adjustedRandIndex(true_labels, cluster_corr(cr, 15))
  }, error = function(e) ari_rmth[i] <<- NA)
}

plot(batch_vals, ari_rmtnet, type = "b", pch = 16, col = "#D62728",
     ylim = c(-0.05, 1), lwd = 2,
     xlab = "Batch effect strength", ylab = "ARI (module recovery)",
     main = "How batch strength kills module recovery")
lines(batch_vals, ari_rmth, type = "b", pch = 17, col = "#1F77B4", lwd = 2)
abline(h = 0, lty = 3, col = "grey50")
legend("topright", c("RMTnet", "Hard threshold"),
       col = c("#D62728", "#1F77B4"), pch = c(16, 17), lwd = 2,
       bty = "n", cex = 0.85)
mtext("RMTnet collapses at batch > 0.2 because the batch eigenvalue\ndominates and the filter can't distinguish it from biology.",
      side = 1, line = 3.5, cex = 0.6, col = "grey30")

# Bar chart: ARI with and without batch correction
scenarios <- c("No batch", "Strong batch\n(no correction)", "Strong batch\n(corrected first)")
rmtnet_aris <- c(
  adjustedRandIndex(true_labels, net_clean$modules),
  adjustedRandIndex(true_labels, net_strong$modules),
  adjustedRandIndex(true_labels,
    rmt_network(mat_corrected, min_module_size = 15, verbose = FALSE)$modules)
)
bp <- barplot(rmtnet_aris, names.arg = scenarios, col = c("#2CA02C", "#D62728", "#1F77B4"),
        border = NA, ylim = c(-0.05, max(rmtnet_aris) + 0.1),
        main = "RMTnet: the fix is simple", ylab = "ARI",
        cex.names = 0.75, las = 1)
text(bp, rmtnet_aris + 0.03, round(rmtnet_aris, 3), cex = 0.8)
mtext("Remove batch effects before running RMTnet (e.g., limma::removeBatchEffect)\nand the modules come back.",
      side = 1, line = 3.5, cex = 0.6, col = "grey30")
dev.off()

cat("All figures saved to", out_dir, "\n\n")

cat("Figure 1: fig1_correlation_heatmaps.png\n")
cat("  Shows what the raw correlation matrix looks like with increasing batch\n\n")
cat("Figure 2: fig2_eigenvalue_spectra.png\n")
cat("  Shows how the batch inflates the largest eigenvalue\n\n")
cat("Figure 3: fig3_spectral_filter_result.png\n")
cat("  Shows the cleaned matrix — batch breaks it, correction fixes it\n\n")
cat("Figure 4: fig4_module_recovery.png\n")
cat("  Shows ARI vs batch strength, and the fix\n")
