# Benchmark Scripts

All scripts used to generate the results reported in the main README. Each can be re-run independently.

**Requirements:** R 4.5+, and the following R packages: `RMThreshold`, `mclust`, `dynamicTreeCut`, `WGCNA`, `msigdbr`, `ALL`, `Biobase`. The probe annotation script requires Python 3 with no external dependencies (uses `urllib` to query MyGene.info).

---

## benchmark.R

**What:** Original clean simulation benchmark — RMTnet vs RMThreshold vs no filtering.

**Data:** `simulate_expression()` with clean block-diagonal structure (10 reps × 5 scenarios varying signal strength, sample size, and number of modules).

**Metric:** Adjusted Rand Index (ARI) against known ground-truth module labels.

**Results:**

| Scenario | RMTnet | RMThreshold | No filtering |
|---|---|---|---|
| Easy (signal = 0.8) | 0.388 | 0.929 | 0.353 |
| Medium (signal = 0.6) | 0.319 | 0.904 | 0.179 |
| Hard (signal = 0.4) | 0.245 | 0.855 | 0.107 |
| Small N (40 samples) | 0.286 | 0.875 | 0.138 |
| Many modules (6) | 0.274 | 0.957 | 0.253 |
| **Overall** | **0.302** | **0.904** | **0.206** |

**Note:** These numbers are from the pre-fix version of the spectral filter. The clean simulation heavily favours hard thresholding because the ground truth is perfectly block-diagonal.

---

## realistic_benchmark.R

**What:** Realistic simulation — overlapping modules, variable gene loadings, unequal module sizes, hub genes, weak signal regimes.

**Data:** Custom `simulate_realistic()` function (10 reps × 5 scenarios). No batch effects — those test preprocessing discipline, not algorithmic quality.

**Metric:** Adjusted Rand Index.

**Results (post spectral filter fix):**

| Scenario | RMTnet | RMThreshold | No filtering |
|---|---|---|---|
| Messy (all features) | 0.243 | 0.931 | 0.183 |
| Heavy overlap (20%) | 0.321 | 0.892 | 0.245 |
| Many hub genes (15) | 0.243 | 0.919 | 0.159 |
| Many small modules | 0.220 | 0.932 | 0.172 |
| Weak signal + noise | 0.092 | 0.748 | 0.074 |
| **Overall** | **0.224** | **0.884** | **0.167** |

**Takeaway:** RMThreshold still wins on ARI because the ground truth is still a discrete partition, even with added complexity.

---

## fuzzy_benchmark.R

**What:** Fuzzy/continuous module simulation — genes load on multiple latent factors with varying weights. No discrete blocks exist. Tests correlation matrix recovery rather than module label recovery.

**Data:** Custom `simulate_fuzzy()` function. Each gene loads continuously on 5–10 latent factors (Dirichlet-like weights, sparsified). 10 reps × 5 scenarios.

**Metrics:**
- Recovery correlation: does the method preserve the ranking of gene-gene correlations? (higher = better)
- Recovery error: Frobenius distance from the true signal matrix, normalised (lower = better)

**Results (post spectral filter fix):**

| Scenario | RMTnet | RMThreshold | Raw |
|---|---|---|---|
| | **Recovery correlation** (higher = better) |
| Standard fuzzy | **0.815** | 0.514 | 0.813 |
| Dense loadings | **0.712** | 0.253 | 0.706 |
| Many factors (10) | **0.736** | 0.417 | 0.731 |
| Few samples | **0.692** | 0.361 | 0.690 |
| Large (1000 genes) | **0.788** | 0.419 | 0.787 |

**Takeaway:** RMTnet preserves the relative ordering of correlations (0.69–0.82) while RMThreshold destroys it (0.25–0.51). Hard thresholding zeroes out real but weak correlations, scrambling which genes are more or less related.

---

## compare.R

**What:** Head-to-head RMTnet vs WGCNA on the same real dataset.

**Data:** ALL leukemia dataset (Chiaretti et al. 2004) — 128 patient samples, 1000 most-variable Affymetrix HGU95Av2 probes, known B-cell vs T-cell subtypes. Installed via Bioconductor `ALL` package.

**Results (post spectral filter fix):**

| Metric | RMTnet | WGCNA |
|---|---|---|
| Modules detected | **26** | 1 |
| Unassigned genes | **3** | 969 |
| Genes assigned | **99.7%** | 3.1% |
| Mean intra-module correlation | 0.311 | 0.678 |
| Best ME \|r\| with B/T label | **0.950** | 0.918 |
| Soft power β | N/A | 14 |
| Runtime | 2.6s | 2.7s |

**Takeaway:** WGCNA collapsed — its R² scale-free topology curve was non-monotonic, leading to an over-aggressive soft power that left 97% of genes unassigned.

---

## enrichment.R

**What:** Pathway enrichment analysis on RMTnet vs WGCNA modules from the ALL dataset.

**Data:** Same ALL leukemia data as `compare.R`. Probe-to-gene-symbol mapping via MyGene.info (see `get_probes.R` + `annotate_probes.py`). Tested against MSigDB Hallmark (50 gene sets) and C7 Immunologic gene sets. Hypergeometric test, BH-corrected FDR < 0.05.

**Results (post spectral filter fix):**

- **RMTnet:** 13 of 26 modules with significant enrichment (50%)
- **WGCNA:** 1 of 1 module with significant enrichment (100% — but it's one generic immune signature)

Selected RMTnet hits: E2F Targets (FDR 7.6×10⁻⁵), MYC Targets (FDR 0.009), TNFa/NF-κB Signalling (FDR 3.3×10⁻¹⁴), KRAS Signalling (FDR 0.027), Interferon Gamma Response (FDR 3.4×10⁻⁵), Heme Metabolism (FDR 2.0×10⁻¹⁰).

**Takeaway:** RMTnet found 13 distinct biological processes. WGCNA found one.

---

## realdata.R

**What:** Standalone RMTnet run on the ALL dataset (no comparison to WGCNA).

**Data:** Same ALL leukemia, 1000 most-variable probes × 128 samples.

**Results:** 26 modules detected, λ₊ = 2.86, 52 signal PCs out of 127 non-zero, Von Neumann entropy 53% of max.

---

## batch_explainer.R

**What:** Visual explanation of why batch effects break spectral filtering. Generates 4 figures showing correlation heatmaps, eigenvalue spectra, spectral filter output, and ARI vs batch strength.

**Data:** Synthetic — 300 genes, 80 samples, 3 modules of 60 genes each, with batch effect strength varied from 0 to 0.5.

**Results:** RMTnet collapses at batch strength > 0.2 because the batch eigenvalue dominates and the filter cannot distinguish it from biology. Hard thresholding is unaffected because the batch raises all correlations uniformly. Fix: run `limma::removeBatchEffect()` before RMTnet — modules come back.

**Output figures:** `fig1_correlation_heatmaps.png`, `fig2_eigenvalue_spectra.png`, `fig3_spectral_filter_result.png`, `fig4_module_recovery.png` (saved to the test directory, not the package).

---

## get_probes.R + annotate_probes.py

**What:** Two-step probe annotation pipeline. `get_probes.R` extracts the top 1000 IQR probe IDs from the ALL dataset. `annotate_probes.py` queries the MyGene.info REST API to map Affymetrix HGU95Av2 probe IDs to HGNC gene symbols.

**Data:** ALL dataset probe IDs → MyGene.info API → gene symbols.

**Results:** 984 of 1000 probes successfully mapped. Output: `probe_symbols.csv` (consumed by `enrichment.R`).

**Why not use hgu95av2.db?** The Bioconductor annotation package `hgu95av2.db` (and its dependency `org.Hs.eg.db`) failed to install on Windows x86 R 4.5.2 due to a persistent SQLite lazy-loading error. The Python API approach is a workaround.
