## compare.R — RMTnet vs WGCNA on ALL leukemia dataset (same 1000 genes, 128 samples)

user_lib <- file.path(Sys.getenv("USERPROFILE"), "R", "library")
.libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(ALL); library(Biobase)
  library(WGCNA)
  library(dynamicTreeCut)
})
disableWGCNAThreads()

pkg_path <- file.path(dirname(getwd()), "test", "RMTnet")
for (f in list.files(file.path(pkg_path, "R"), full.names = TRUE)) source(f)

data(ALL)

# ── Shared data prep (identical for both methods) ────────────────────────────
mat_full <- exprs(ALL)
iqr      <- apply(mat_full, 1, IQR)
top1000  <- order(iqr, decreasing = TRUE)[1:1000]
mat      <- mat_full[top1000, ]          # 1000 genes x 128 samples
datExpr  <- t(mat)                       # WGCNA wants samples x genes

cell_type <- substr(as.character(ALL$BT), 1, 1)  # "B" or "T"
bt_numeric <- as.integer(cell_type == "T")        # 0=B, 1=T

# ── Helper: module eigengene correlation with B/T label ──────────────────────
best_me_cor <- function(expr_mat, modules) {
  # expr_mat: genes x samples; modules: integer vector per gene (0 = unassigned)
  mods <- sort(unique(modules[modules > 0]))
  if (length(mods) == 0) return(NA_real_)
  cors <- sapply(mods, function(m) {
    g    <- which(modules == m)
    if (length(g) < 2) return(NA_real_)
    me   <- prcomp(t(expr_mat[g, ]), rank. = 1)$x[, 1]
    abs(cor(me, bt_numeric, use = "complete.obs"))
  })
  max(cors, na.rm = TRUE)
}

# ── Helper: mean intra-module correlation ────────────────────────────────────
mean_intra_cor <- function(expr_mat, modules, corr_mat) {
  mods <- sort(unique(modules[modules > 0]))
  if (length(mods) == 0) return(NA_real_)
  vals <- sapply(mods, function(m) {
    g <- which(modules == m)
    if (length(g) < 2) return(NA_real_)
    sub <- corr_mat[g, g]
    mean(sub[upper.tri(sub)])
  })
  mean(vals, na.rm = TRUE)
}

# ═══════════════════════════════════════════════════════════════════════════
# METHOD 1: RMTnet
# ═══════════════════════════════════════════════════════════════════════════
cat("Running RMTnet...\n")
t0_rmt <- proc.time()["elapsed"]
net    <- rmt_network(mat, min_module_size = 15, verbose = FALSE)
t1_rmt <- proc.time()["elapsed"]

rmt_modules    <- net$modules                         # integer, 0 = unassigned
rmt_n_mod      <- max(rmt_modules)
rmt_unassigned <- sum(rmt_modules == 0)
rmt_intra_cor  <- mean_intra_cor(mat, rmt_modules, net$cleaned_corr)
rmt_me_cor     <- best_me_cor(mat, rmt_modules)
rmt_time       <- round(t1_rmt - t0_rmt, 1)

cat(sprintf("  Done: %d modules, %d unassigned genes (%.1fs)\n\n",
            rmt_n_mod, rmt_unassigned, rmt_time))

# ═══════════════════════════════════════════════════════════════════════════
# METHOD 2: WGCNA
# ═══════════════════════════════════════════════════════════════════════════
cat("Running WGCNA...\n")
cat("  Step 1: picking soft threshold (testing powers 1-20)...\n")
t0_wgcna <- proc.time()["elapsed"]

sft <- pickSoftThreshold(datExpr, powerVector = 1:20,
                         networkType = "signed hybrid",
                         verbose = 0)

# Pick lowest power where R² >= 0.85, fallback to best available
r2  <- sft$fitIndices[, 2]
ok  <- which(r2 >= 0.85)
soft_power <- if (length(ok) > 0) sft$fitIndices[min(ok), 1] else sft$fitIndices[which.max(r2), 1]
cat(sprintf("  Soft power chosen: %d  (R² = %.3f)\n", soft_power, max(r2, na.rm=TRUE)))

cat("  Step 2: blockwiseModules...\n")
bwm <- blockwiseModules(
  datExpr,
  power            = soft_power,
  networkType      = "signed hybrid",
  minModuleSize    = 15,
  mergeCutHeight   = 0.25,
  numericLabels    = TRUE,
  verbose          = 0
)
t1_wgcna <- proc.time()["elapsed"]

wgcna_colors   <- bwm$colors                          # 0 = grey (unassigned)
wgcna_n_mod    <- length(unique(wgcna_colors[wgcna_colors > 0]))
wgcna_unassign <- sum(wgcna_colors == 0)
raw_corr       <- cor(t(mat))
wgcna_intra    <- mean_intra_cor(mat, wgcna_colors, raw_corr)
wgcna_me_cor   <- best_me_cor(mat, wgcna_colors)
wgcna_time     <- round(t1_wgcna - t0_wgcna, 1)

cat(sprintf("  Done: %d modules, %d unassigned genes (%.1fs)\n\n",
            wgcna_n_mod, wgcna_unassign, wgcna_time))

# ═══════════════════════════════════════════════════════════════════════════
# COMPARISON TABLE
# ═══════════════════════════════════════════════════════════════════════════
cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║          RMTnet vs WGCNA — ALL Leukemia Dataset          ║\n")
cat("║          1000 genes  |  128 samples  |  B/T-cell         ║\n")
cat("╠════════════════════════════════════╦══════════╦══════════╣\n")
cat("║ Metric                             ║  RMTnet  ║  WGCNA   ║\n")
cat("╠════════════════════════════════════╬══════════╬══════════╣\n")
cat(sprintf("║ Modules detected                   ║  %6d  ║  %6d  ║\n", rmt_n_mod,      wgcna_n_mod))
cat(sprintf("║ Unassigned genes                   ║  %6d  ║  %6d  ║\n", rmt_unassigned,  wgcna_unassign))
cat(sprintf("║ Pct assigned                       ║  %5.1f%%  ║  %5.1f%%  ║\n",
            100*(1000-rmt_unassigned)/1000, 100*(1000-wgcna_unassign)/1000))
cat(sprintf("║ Mean intra-module correlation      ║  %6.3f  ║  %6.3f  ║\n", rmt_intra_cor,  wgcna_intra))
cat(sprintf("║ Best module eigengene |r| w/ B/T   ║  %6.3f  ║  %6.3f  ║\n", rmt_me_cor,     wgcna_me_cor))
cat(sprintf("║ Soft power used                    ║    N/A   ║  %6d  ║\n", soft_power))
cat(sprintf("║ Runtime (sec)                      ║  %6.1f  ║  %6.1f  ║\n", rmt_time,       wgcna_time))
cat("╚════════════════════════════════════╩══════════╩══════════╝\n")

cat("\nNotes:\n")
cat("  Mean intra-module correlation: average Pearson r of gene pairs within same module\n")
cat("  Best ME |r| with B/T: max |correlation| between any module eigengene and B/T label\n")
cat("  RMTnet uses cleaned (noise-filtered) correlation; WGCNA uses raw correlation\n")
