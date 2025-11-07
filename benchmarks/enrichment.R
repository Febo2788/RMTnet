## enrichment.R — pathway enrichment on RMTnet vs WGCNA modules (ALL leukemia)

out_file <- "C:/Users/felix/downloads/test/enrichment_results.txt"
sink(out_file, split = TRUE)  # split=TRUE means output goes to both console and file

user_lib <- file.path(Sys.getenv("USERPROFILE"), "R", "library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

# Install what we need
bioc_pkgs <- c("ALL", "Biobase")
cran_pkgs <- c("msigdbr", "dynamicTreeCut")
for (p in bioc_pkgs)
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p, lib = user_lib, ask = FALSE, update = FALSE)
for (p in cran_pkgs)
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, lib = user_lib, repos = "https://cloud.r-project.org", quiet = TRUE)

suppressPackageStartupMessages({
  library(ALL); library(Biobase)
  library(msigdbr)
  library(dynamicTreeCut)
  library(WGCNA)
})
disableWGCNAThreads()

pkg_path <- file.path(dirname(getwd()), "test", "RMTnet")
for (f in list.files(file.path(pkg_path, "R"), full.names = TRUE)) source(f)

# ── Data prep (identical to previous runs) ───────────────────────────────────
data(ALL)
mat_full <- exprs(ALL)
iqr      <- apply(mat_full, 1, IQR)
top_idx  <- order(iqr, decreasing = TRUE)[1:1000]
mat      <- mat_full[top_idx, ]

# ── Probe → gene symbol mapping (from Python/MyGene.info lookup) ─────────────
annot     <- read.csv("C:/Users/felix/downloads/test/probe_symbols.csv",
                      stringsAsFactors = FALSE)
sym_map   <- setNames(annot$symbol, annot$probe_id)
sym_map   <- sym_map[sym_map != ""]
cat(sprintf("Probes mapped to gene symbols: %d / %d\n\n",
            length(sym_map), nrow(annot)))

# ── MSigDB gene sets ──────────────────────────────────────────────────────────
# Hallmark (H): 50 well-curated biological states
# C7 immunologic: relevant for leukemia
h_sets  <- msigdbr(species = "Homo sapiens", category = "H")
c7_sets <- msigdbr(species = "Homo sapiens", category = "C7",
                   subcategory = "IMMUNESIGDB")

make_gset_list <- function(df) {
  split(df$gene_symbol, df$gs_name)
}
gsets_h  <- make_gset_list(h_sets)
gsets_c7 <- make_gset_list(c7_sets)

# ── Hypergeometric enrichment ─────────────────────────────────────────────────
# universe = all mapped probes in our 1000-gene set
universe <- unique(sym_map)
N        <- length(universe)

enrich_module <- function(module_genes, gset_list, universe, fdr_cutoff = 0.05) {
  mod_syms <- unique(module_genes[!is.na(module_genes)])
  if (length(mod_syms) < 3) return(NULL)

  res <- lapply(names(gset_list), function(gs_name) {
    gs      <- intersect(gset_list[[gs_name]], universe)
    overlap <- length(intersect(mod_syms, gs))
    if (overlap == 0) return(NULL)
    p <- phyper(overlap - 1, length(gs), N - length(gs),
                length(mod_syms), lower.tail = FALSE)
    data.frame(pathway = gs_name, overlap = overlap,
               mod_size = length(mod_syms), gs_size = length(gs),
               p_val = p, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, Filter(Negate(is.null), res))
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res$fdr <- p.adjust(res$p_val, method = "BH")
  res <- res[order(res$fdr), ]
  res[res$fdr < fdr_cutoff, , drop = FALSE]
}

# ── Run RMTnet ────────────────────────────────────────────────────────────────
cat("Running RMTnet...\n")
net        <- rmt_network(mat, min_module_size = 15, verbose = FALSE)
rmt_mods   <- net$modules

# ── Run WGCNA ────────────────────────────────────────────────────────────────
cat("Running WGCNA...\n")
sft        <- pickSoftThreshold(t(mat), powerVector = 1:20,
                                networkType = "signed hybrid", verbose = 0)
r2         <- sft$fitIndices[, 2]
soft_power <- sft$fitIndices[which(r2 >= 0.85)[1], 1]
bwm        <- blockwiseModules(t(mat), power = soft_power,
                               networkType = "signed hybrid",
                               minModuleSize = 15, mergeCutHeight = 0.25,
                               numericLabels = TRUE, verbose = 0)
wgcna_mods <- bwm$colors

# ── Enrichment per module ─────────────────────────────────────────────────────

run_enrichment <- function(modules, sym_map, gsets_h, gsets_c7, label) {
  mod_ids  <- sort(unique(modules[modules > 0]))
  cat(sprintf("\n%s — %d modules\n", label, length(mod_ids)))
  cat(strrep("─", 70), "\n")

  all_hits <- list()
  for (m in mod_ids) {
    probe_in_mod  <- names(modules)[modules == m]
    genes_in_mod  <- sym_map[probe_in_mod]

    hits_h  <- enrich_module(genes_in_mod, gsets_h,  universe)
    hits_c7 <- enrich_module(genes_in_mod, gsets_c7, universe)

    top_h  <- if (!is.null(hits_h)  && nrow(hits_h)  > 0) hits_h[1, ]  else NULL
    top_c7 <- if (!is.null(hits_c7) && nrow(hits_c7) > 0) hits_c7[1, ] else NULL

    if (!is.null(top_h) || !is.null(top_c7)) {
      cat(sprintf("\nModule %2d  (%d genes)\n", m, sum(modules == m)))
      if (!is.null(top_h))
        cat(sprintf("  Hallmark: %-55s  overlap=%d  FDR=%.3g\n",
                    gsub("HALLMARK_", "", top_h$pathway),
                    top_h$overlap, top_h$fdr))
      if (!is.null(top_c7))
        cat(sprintf("  Immunol:  %-55s  overlap=%d  FDR=%.3g\n",
                    substr(top_c7$pathway, 1, 55),
                    top_c7$overlap, top_c7$fdr))
      all_hits[[as.character(m)]] <- list(hallmark = top_h, c7 = top_c7)
    }
  }

  n_enriched <- length(all_hits)
  cat(sprintf("\n%s: %d / %d modules with significant enrichment (FDR < 0.05)\n",
              label, n_enriched, length(mod_ids)))
  invisible(all_hits)
}

rmt_hits   <- run_enrichment(rmt_mods,   sym_map, gsets_h, gsets_c7, "RMTnet")
wgcna_hits <- run_enrichment(wgcna_mods, sym_map, gsets_h, gsets_c7, "WGCNA")

# ── Summary table ─────────────────────────────────────────────────────────────
cat("\n")
cat("╔═══════════════════════════════════════════════════╗\n")
cat("║       Enrichment Summary — ALL Leukemia           ║\n")
cat("╠══════════════════════════════╦════════╦═══════════╣\n")
cat("║ Metric                       ║ RMTnet ║   WGCNA   ║\n")
cat("╠══════════════════════════════╬════════╬═══════════╣\n")
cat(sprintf("║ Modules detected             ║  %4d  ║  %7d  ║\n",
            max(rmt_mods), max(wgcna_mods)))
cat(sprintf("║ Modules with FDR<0.05 hit    ║  %4d  ║  %7d  ║\n",
            length(rmt_hits), length(wgcna_hits)))
cat(sprintf("║ Enrichment rate              ║  %3.0f%%  ║  %7s  ║\n",
            100 * length(rmt_hits) / max(rmt_mods),
            paste0(round(100 * length(wgcna_hits) / max(max(wgcna_mods),1)), "%")))
cat("╚══════════════════════════════╩════════╩═══════════╝\n")
sink()
cat("\nResults saved to", out_file, "\n")
