## realdata.R — run RMTnet on the ALL leukemia dataset
## ALL: 128 patients, ~12,600 Affymetrix probes, known B/T-cell subtypes

user_lib <- file.path(Sys.getenv("USERPROFILE"), "R", "library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

# Install ALL from Bioconductor if needed
if (!requireNamespace("ALL", quietly = TRUE)) {
  BiocManager::install("ALL", lib = user_lib, ask = FALSE, quiet = TRUE)
}

library(ALL)
data(ALL)

# ── Load RMTnet ──────────────────────────────────────────────────────────────
pkg_path <- file.path(dirname(getwd()), "test", "RMTnet")
for (f in list.files(file.path(pkg_path, "R"), full.names = TRUE)) source(f)
library(dynamicTreeCut)

# ── Prep expression matrix ───────────────────────────────────────────────────
cat("Dataset: ALL leukemia (Chiaretti et al. 2004)\n")
cat("Samples:", ncol(ALL), "\n")
cat("Probes: ", nrow(ALL), "\n\n")

# Take top 1000 most variable probes (IQR-based, standard practice)
mat_full <- exprs(ALL)
iqr      <- apply(mat_full, 1, IQR)
top1000  <- order(iqr, decreasing = TRUE)[1:1000]
mat      <- mat_full[top1000, ]

cat("Using top 1000 most variable probes\n")
cat("Matrix: ", nrow(mat), "genes x", ncol(mat), "samples\n")
cat("Q = T/N =", round(ncol(mat)/nrow(mat), 3), "\n\n")

# ── Run RMTnet ───────────────────────────────────────────────────────────────
cat("Running RMTnet...\n")
net <- rmt_network(mat, min_module_size = 15, verbose = TRUE)

cat("\n")
print(net)

# ── Compare modules to known biology ─────────────────────────────────────────
cat("\n══════════════════════════════════════════\n")
cat("Comparing modules to known phenotypes\n")
cat("══════════════════════════════════════════\n\n")

# B-cell vs T-cell (BT column: B, B1, B2, B3, B4, T, T1, T2, T3, T4)
cell_type  <- substr(as.character(ALL$BT), 1, 1)  # "B" or "T"
mol_remiss <- as.character(ALL$mol.biol)           # BCR/ABL, ALL1/AF4, NEG, etc.

cat("Cell type distribution:\n")
print(table(cell_type))
cat("\nMolecular biology distribution:\n")
print(table(mol_remiss))

# For each module, show enrichment for B vs T
cat("\n--- Module x Cell type table ---\n")
mod_labels <- net$modules
named_mods <- mod_labels[mod_labels > 0]  # exclude unassigned (0)
named_ct   <- cell_type[mod_labels > 0]
if (length(unique(named_mods)) > 0) {
  print(table(Module = named_mods, CellType = named_ct))
} else {
  cat("No assigned modules\n")
}

cat("\n--- Module table ---\n")
print(net$module_table)

cat("\n--- Key RMT statistics ---\n")
cat(sprintf("Signal PCs (lambda > lambda+): %d / %d (%.1f%%)\n",
            net$n_signal, nrow(mat), 100 * net$n_signal / nrow(mat)))
cat(sprintf("MP noise ceiling lambda+:      %.4f\n", net$lambda_plus))
cat(sprintf("Von Neumann entropy:           %.3f / %.3f (%.1f%% of max)\n",
            net$entropy, net$entropy_max,
            100 * net$entropy / net$entropy_max))
