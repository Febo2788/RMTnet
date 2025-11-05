#' Clean a gene correlation matrix by removing noise eigenvalues
#'
#' Applies spectral filtering based on the Marchenko-Pastur threshold.
#' Noise eigenvalues (those below \code{lambda_plus}) are replaced with their
#' mean, preserving the matrix trace (= N). Signal eigenvalues are kept
#' unchanged. The cleaned matrix is then reconstructed via eigen-decomposition.
#'
#' This is analogous to \code{limma::removeBatchEffect()} but instead of
#' removing known confounders it removes statistically proven noise — the
#' part of your correlation matrix that would exist even if all your genes
#' were completely independent, purely due to finite sample size.
#'
#' @param mp_result The output of \code{\link{mp_threshold}}.
#'
#' @return A list with:
#'   \describe{
#'     \item{cleaned_corr}{The noise-filtered correlation matrix (N x N).}
#'     \item{raw_corr}{The original unfiltered correlation matrix.}
#'     \item{lambda_signal}{Signal eigenvalues (above threshold).}
#'     \item{lambda_noise_mean}{Mean of removed noise eigenvalues.}
#'     \item{n_signal}{Number of retained signal eigenvalues.}
#'     \item{mp_threshold}{The lambda_plus cutoff used.}
#'   }
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rnorm(200 * 50), nrow = 200, ncol = 50)
#' rownames(mat) <- paste0("gene", seq_len(200))
#' mp  <- mp_threshold(mat)
#' filtered <- spectral_filter(mp)
#' dim(filtered$cleaned_corr)
#'
#' @export
spectral_filter <- function(mp_result) {

  .validate_mp_result(mp_result)

  eigenvalues  <- mp_result$eigenvalues
  eigenvectors <- mp_result$eigenvectors
  n_signal     <- mp_result$n_signal
  lambda_plus  <- mp_result$lambda_plus
  N            <- mp_result$N

  # Replace noise eigenvalues with their mean (preserves trace = N)
  noise_mean   <- if (n_signal < N) mean(eigenvalues[(n_signal + 1L):N]) else 0
  lam_clean    <- eigenvalues
  if (n_signal < N) lam_clean[(n_signal + 1L):N] <- noise_mean

  # Reconstruct:  C_clean = U * diag(lam_clean) * U^T
  C_clean <- eigenvectors %*% diag(lam_clean) %*% t(eigenvectors)

  # Re-normalise diagonal to exactly 1 (numerical precision)
  d       <- sqrt(diag(C_clean))
  C_clean <- C_clean / outer(d, d)
  diag(C_clean) <- 1

  # Carry over gene names if present
  genes <- rownames(mp_result$corr_matrix)
  if (!is.null(genes)) {
    rownames(C_clean) <- colnames(C_clean) <- genes
    rownames(mp_result$corr_matrix) <- colnames(mp_result$corr_matrix) <- genes
  }

  list(
    cleaned_corr     = C_clean,
    raw_corr         = mp_result$corr_matrix,
    lambda_signal    = eigenvalues[seq_len(n_signal)],
    lambda_noise_mean= noise_mean,
    n_signal         = n_signal,
    mp_threshold     = lambda_plus
  )
}


#' @noRd
.validate_mp_result <- function(x) {
  needed <- c("eigenvalues", "eigenvectors", "n_signal",
              "lambda_plus", "N", "corr_matrix")
  missing <- setdiff(needed, names(x))
  if (length(missing))
    stop("mp_result is missing fields: ", paste(missing, collapse = ", "),
         ". Did you run mp_threshold() first?")
}
