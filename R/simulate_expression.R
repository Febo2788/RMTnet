#' Simulate a gene expression matrix with embedded co-expression modules
#'
#' Generates a synthetic expression matrix useful for testing and vignettes.
#' Genes within the same module share a latent factor, producing genuine
#' co-expression structure that RMT should detect above the noise floor.
#'
#' @param n_genes Total number of genes. Default \code{500}.
#' @param n_samples Number of samples. Default \code{80}.
#' @param n_modules Number of embedded co-expression modules. Default \code{3}.
#' @param module_size Number of genes per module. Default \code{50}.
#' @param signal_strength Strength of the within-module signal (0-1).
#'   Higher values = stronger co-expression. Default \code{0.7}.
#' @param seed Random seed for reproducibility. Default \code{42}.
#'
#' @return A numeric matrix (n_genes x n_samples) with rownames
#'   \code{"gene1"}, \code{"gene2"}, etc.
#'
#' @examples
#' mat <- simulate_expression(n_genes = 300, n_samples = 60, n_modules = 3)
#' dim(mat)
#'
#' @export
simulate_expression <- function(n_genes         = 500L,
                                 n_samples        = 80L,
                                 n_modules        = 3L,
                                 module_size      = 50L,
                                 signal_strength  = 0.7,
                                 seed             = 42L) {
  set.seed(seed)

  mat <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)

  # Inject latent factors for each module
  for (m in seq_len(n_modules)) {
    start <- (m - 1L) * module_size + 1L
    end   <- min(m * module_size, n_genes)
    if (start > n_genes) break

    latent <- rnorm(n_samples)   # shared signal for this module
    for (g in seq(start, end)) {
      mat[g, ] <- sqrt(1 - signal_strength) * mat[g, ] +
                  sqrt(signal_strength)     * latent
    }
  }

  rownames(mat) <- paste0("gene", seq_len(n_genes))
  colnames(mat) <- paste0("sample", seq_len(n_samples))
  mat
}
