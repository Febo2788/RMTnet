#' Fit the Marchenko-Pastur noise threshold for a matrix
#'
#' Given an expression or correlation matrix, computes the Marchenko-Pastur
#' upper bound (lambda_plus). Eigenvalues of the correlation matrix that exceed
#' this threshold are considered statistically significant signal; those below
#' are indistinguishable from random noise given your sample size.
#'
#' The Marchenko-Pastur law describes the limiting spectral distribution of a
#' Wishart random matrix W = (1/T) X Xt, where X has i.i.d. entries. The
#' bulk eigenvalue support is bounded by:
#'
#'   lambda_plus  = sigma^2 * (1 + 1/sqrt(Q))^2
#'   lambda_minus = sigma^2 * (1 - 1/sqrt(Q))^2
#'
#' where Q = T/N (samples / genes) and sigma^2 is fitted by maximum likelihood.
#'
#' @param expr_matrix A numeric matrix of shape (genes x samples), e.g.
#'   log-normalised counts from DESeq2 or edgeR. Rows are features, columns
#'   are samples. Can also accept a SummarizedExperiment object.
#' @param assay_name If a SummarizedExperiment is supplied, the assay to use.
#'   Default \code{"counts"}.
#'
#' @return A list with:
#'   \describe{
#'     \item{lambda_plus}{The MP noise ceiling. Eigenvalues above this are signal.}
#'     \item{lambda_minus}{The MP noise floor.}
#'     \item{sigma2}{Fitted noise variance.}
#'     \item{Q}{The matrix ratio T/N.}
#'     \item{N}{Number of features (genes).}
#'     \item{T}{Number of observations (samples).}
#'     \item{eigenvalues}{All eigenvalues of the correlation matrix, descending.}
#'     \item{n_signal}{Number of eigenvalues above lambda_plus.}
#'     \item{n_noise}{Number of eigenvalues below lambda_plus.}
#'   }
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rnorm(200 * 50), nrow = 200, ncol = 50)
#' rownames(mat) <- paste0("gene", seq_len(200))
#' result <- mp_threshold(mat)
#' result$n_signal
#' result$lambda_plus
#'
#' @export
mp_threshold <- function(expr_matrix, assay_name = "counts") {

  expr_matrix <- .extract_matrix(expr_matrix, assay_name)
  .validate_matrix(expr_matrix)

  N <- nrow(expr_matrix)   # genes
  T <- ncol(expr_matrix)   # samples
  Q <- T / N

  # Standardise each gene (zero mean, unit variance) before correlating
  expr_std <- t(scale(t(expr_matrix)))

  # Pearson correlation matrix  C = (1/T) X Xt
  C <- tcrossprod(expr_std) / T

  # Eigendecomposition (symmetric matrix — use efficient eigen())
  eig      <- eigen(C, symmetric = TRUE)
  eigenvalues <- eig$values   # already descending

  # Fit sigma^2 by maximising the MP log-likelihood over the bulk
  sigma2 <- .fit_mp_sigma2(eigenvalues, Q)

  lambda_plus  <- sigma2 * (1 + 1 / sqrt(Q))^2
  lambda_minus <- sigma2 * (1 - 1 / sqrt(Q))^2

  n_signal <- sum(eigenvalues > lambda_plus)
  n_noise  <- N - n_signal

  list(
    lambda_plus  = lambda_plus,
    lambda_minus = lambda_minus,
    sigma2       = sigma2,
    Q            = Q,
    N            = N,
    T            = T,
    eigenvalues  = eigenvalues,
    eigenvectors = eig$vectors,
    corr_matrix  = C,
    n_signal     = n_signal,
    n_noise      = n_noise
  )
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Marchenko-Pastur density at given eigenvalues
#' @noRd
.mp_pdf <- function(lambda, sigma2, Q) {
  lp <- sigma2 * (1 + 1 / sqrt(Q))^2
  lm <- sigma2 * (1 - 1 / sqrt(Q))^2
  in_bulk <- lambda >= lm & lambda <= lp
  rho <- numeric(length(lambda))
  rho[in_bulk] <- (Q / (2 * pi * sigma2)) *
    sqrt(pmax((lp - lambda[in_bulk]) * (lambda[in_bulk] - lm), 0)) /
    lambda[in_bulk]
  rho
}

#' Fit sigma^2 by maximising the MP log-likelihood
#' @noRd
.fit_mp_sigma2 <- function(eigenvalues, Q) {
  # Search over log(sigma^2) for numerical stability
  neg_ll <- function(log_s2) {
    s2 <- exp(log_s2)
    lp <- s2 * (1 + 1 / sqrt(Q))^2
    bulk <- eigenvalues[eigenvalues <= lp]
    if (length(bulk) < 5L) return(1e9)
    pdf_vals <- .mp_pdf(bulk, s2, Q)
    pdf_vals <- pmax(pdf_vals, 1e-15)
    -sum(log(pdf_vals))
  }
  result <- optimise(neg_ll, interval = c(-3, 3))
  exp(result$minimum)
}

#' Extract a plain numeric matrix from various input types
#' @noRd
.extract_matrix <- function(x, assay_name) {
  if (inherits(x, "SummarizedExperiment")) {
    return(as.matrix(SummarizedExperiment::assay(x, assay_name)))
  }
  if (inherits(x, "matrix")) return(x)
  if (inherits(x, "data.frame")) return(as.matrix(x))
  stop("expr_matrix must be a matrix, data.frame, or SummarizedExperiment.")
}

#' Basic sanity checks on the input matrix
#' @noRd
.validate_matrix <- function(mat) {
  if (!is.numeric(mat))
    stop("expr_matrix must be numeric.")
  if (nrow(mat) < 10L)
    stop("Need at least 10 genes for meaningful RMT analysis.")
  if (ncol(mat) < 10L)
    stop("Need at least 10 samples for meaningful RMT analysis.")
  if (ncol(mat) >= nrow(mat))
    warning(
      "More samples than genes (T >= N). RMT assumes N >> T is not required, ",
      "but results are most interpretable when N > T."
    )
  if (any(!is.finite(mat)))
    stop("expr_matrix contains NA, NaN, or Inf values. Please impute or filter first.")
  invisible(NULL)
}
