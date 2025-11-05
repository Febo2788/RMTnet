#' Build a gene co-expression network using Random Matrix Theory
#'
#' The main entry point for RMTnet. Takes an expression matrix and returns a
#' cleaned co-expression network with gene modules — a parameter-free
#' alternative to WGCNA.
#'
#' ## How it differs from WGCNA
#'
#' WGCNA requires you to manually select a soft-thresholding power beta by
#' inspecting a scale-free topology plot. The choice is subjective and changes
#' your network substantially. RMTnet derives the noise threshold mathematically
#' from your data dimensions using the Marchenko-Pastur law — no parameters
#' to tune, reproducible results every time.
#'
#' ## Pipeline
#' \enumerate{
#'   \item Standardise expression matrix (zero mean, unit variance per gene)
#'   \item Compute Pearson correlation matrix
#'   \item Fit Marchenko-Pastur distribution to find noise ceiling (lambda_plus)
#'   \item Remove noise eigenvalues via spectral filtering
#'   \item Detect modules via hierarchical clustering on cleaned correlations
#' }
#'
#' @param expr_matrix A numeric matrix (genes x samples) or SummarizedExperiment.
#'   Should be normalised (e.g. VST, rlog, or log-CPM). Rows = genes,
#'   columns = samples.
#' @param assay_name Assay to use if a SummarizedExperiment is supplied.
#'   Default \code{"counts"}.
#' @param min_module_size Minimum number of genes per module. Default \code{30}.
#'   Passed to \code{dynamicTreeCut::cutreeDynamic}.
#' @param cor_threshold After cleaning, edges below this absolute correlation
#'   are set to zero in the adjacency matrix. Default \code{0.3}. Set to
#'   \code{0} to keep all edges.
#' @param verbose Print progress messages. Default \code{TRUE}.
#'
#' @return An object of class \code{"RMTnetwork"}, which is a list with:
#'   \describe{
#'     \item{modules}{Named integer vector: gene -> module number (0 = unassigned).}
#'     \item{module_table}{data.frame summarising each module (size, mean correlation).}
#'     \item{cleaned_corr}{The noise-filtered N x N correlation matrix.}
#'     \item{adjacency}{Thresholded adjacency matrix used for module detection.}
#'     \item{mp}{Raw output of \code{mp_threshold()} for diagnostics.}
#'     \item{filter}{Raw output of \code{spectral_filter()} for diagnostics.}
#'     \item{n_signal}{Number of statistically significant PCs.}
#'     \item{lambda_plus}{The Marchenko-Pastur noise ceiling.}
#'     \item{entropy}{Von Neumann entropy of the correlation matrix (nats).}
#'     \item{entropy_max}{Maximum possible entropy log(N).}
#'     \item{call}{The matched call.}
#'   }
#'
#' @examples
#' set.seed(42)
#' # Simulate 500 genes, 80 samples with 3 embedded modules
#' mat <- simulate_expression(n_genes = 500, n_samples = 80, n_modules = 3)
#' net <- rmt_network(mat)
#' print(net)
#' plot(net)
#'
#' @seealso \code{\link{mp_threshold}}, \code{\link{spectral_filter}},
#'   \code{\link{plot.RMTnetwork}}, \code{\link{simulate_expression}}
#'
#' @export
rmt_network <- function(expr_matrix,
                         assay_name     = "counts",
                         min_module_size = 30L,
                         cor_threshold  = 0.3,
                         verbose        = TRUE) {

  call <- match.call()

  .msg(verbose, "Step 1/4  Computing correlation matrix and MP threshold...")
  mp <- mp_threshold(expr_matrix, assay_name = assay_name)

  .msg(verbose, sprintf(
    "         N=%d genes, T=%d samples, Q=%.2f",
    mp$N, mp$T, mp$Q
  ))
  .msg(verbose, sprintf(
    "         lambda_plus=%.4f  |  signal PCs: %d / %d  (%.0f%%)",
    mp$lambda_plus, mp$n_signal, mp$N, 100 * mp$n_signal / mp$N
  ))

  .msg(verbose, "Step 2/4  Removing noise eigenvalues (spectral filtering)...")
  filt <- spectral_filter(mp)

  .msg(verbose, "Step 3/4  Building adjacency and detecting modules...")
  adj <- filt$cleaned_corr
  adj[abs(adj) < cor_threshold] <- 0
  diag(adj) <- 0

  # Convert correlations to distances for clustering
  dist_mat <- as.dist(1 - abs(filt$cleaned_corr))
  tree     <- hclust(dist_mat, method = "average")

  if (!requireNamespace("dynamicTreeCut", quietly = TRUE)) {
    stop("Package 'dynamicTreeCut' is required. Install with:\n",
         "  install.packages('dynamicTreeCut')")
  }

  modules <- dynamicTreeCut::cutreeDynamic(
    dendro          = tree,
    distM           = as.matrix(dist_mat),
    deepSplit       = 2,
    minClusterSize  = min_module_size,
    verbose         = 0L
  )
  names(modules) <- rownames(filt$cleaned_corr)

  .msg(verbose, "Step 4/4  Computing summary statistics...")
  entropy     <- .von_neumann_entropy(mp$eigenvalues)
  entropy_max <- log(mp$N)

  module_table <- .summarise_modules(modules, filt$cleaned_corr)

  .msg(verbose, sprintf(
    "Done. %d modules detected  |  Von Neumann entropy: %.3f / %.3f (%.0f%%)",
    nrow(module_table), entropy, entropy_max, 100 * entropy / entropy_max
  ))

  structure(
    list(
      modules      = modules,
      module_table = module_table,
      cleaned_corr = filt$cleaned_corr,
      adjacency    = adj,
      mp           = mp,
      filter       = filt,
      n_signal     = mp$n_signal,
      lambda_plus  = mp$lambda_plus,
      entropy      = entropy,
      entropy_max  = entropy_max,
      call         = call
    ),
    class = "RMTnetwork"
  )
}


# ── S3 methods ────────────────────────────────────────────────────────────────

#' Print an RMTnetwork object
#' @param x An RMTnetwork object.
#' @param ... Ignored.
#' @export
print.RMTnetwork <- function(x, ...) {
  cat("── RMTnet Co-Expression Network ──────────────────────────\n")
  cat(sprintf("  Genes:            %d\n", x$mp$N))
  cat(sprintf("  Samples:          %d\n", x$mp$T))
  cat(sprintf("  Matrix ratio Q:   %.2f\n", x$mp$Q))
  cat(sprintf("  MP threshold λ+:  %.4f\n", x$lambda_plus))
  cat(sprintf("  Signal PCs:       %d / %d  (%.0f%% real)\n",
              x$n_signal, x$mp$N, 100 * x$n_signal / x$mp$N))
  cat(sprintf("  Von Neumann S:    %.3f / %.3f nats  (%.0f%% of max)\n",
              x$entropy, x$entropy_max, 100 * x$entropy / x$entropy_max))
  cat(sprintf("  Modules detected: %d\n", nrow(x$module_table)))
  cat("──────────────────────────────────────────────────────────\n")
  if (nrow(x$module_table) > 0) {
    cat("  Module sizes:\n")
    print(x$module_table, row.names = FALSE)
  }
  invisible(x)
}

#' Summary of an RMTnetwork object
#' @param object An RMTnetwork object.
#' @param ... Ignored.
#' @export
summary.RMTnetwork <- function(object, ...) {
  print(object)
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Von Neumann entropy of the eigenvalue spectrum
#' @noRd
.von_neumann_entropy <- function(eigenvalues) {
  p  <- eigenvalues / sum(eigenvalues)
  p  <- p[p > 0]
  -sum(p * log(p))
}

#' Build module summary table
#' @noRd
.summarise_modules <- function(modules, corr_matrix) {
  mod_ids <- sort(unique(modules[modules != 0]))
  if (length(mod_ids) == 0L) {
    return(data.frame(module = integer(0), size = integer(0),
                      mean_cor = numeric(0)))
  }
  rows <- lapply(mod_ids, function(m) {
    genes   <- names(modules)[modules == m]
    sub_cor <- corr_matrix[genes, genes, drop = FALSE]
    # Mean off-diagonal correlation
    off_diag <- sub_cor[upper.tri(sub_cor)]
    data.frame(
      module   = m,
      size     = length(genes),
      mean_cor = round(mean(off_diag), 4)
    )
  })
  do.call(rbind, rows)
}

#' Conditional message printer
#' @noRd
.msg <- function(verbose, ...) {
  if (verbose) message(...)
}
