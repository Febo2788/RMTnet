#' Plot diagnostics for an RMTnetwork object
#'
#' Produces a multi-panel diagnostic figure:
#' \enumerate{
#'   \item Eigenvalue spectrum with Marchenko-Pastur bounds
#'   \item Marchenko-Pastur PDF fit to empirical bulk eigenvalues
#'   \item Raw vs cleaned correlation heatmaps
#'   \item Module size distribution
#' }
#'
#' @param x An \code{RMTnetwork} object from \code{\link{rmt_network}}.
#' @param which Which panels to show. One of \code{"all"}, \code{"spectrum"},
#'   \code{"mp_fit"}, \code{"heatmap"}, or \code{"modules"}. Default \code{"all"}.
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}. Called for side effects (plots).
#'
#' @examples
#' set.seed(42)
#' mat <- simulate_expression(n_genes = 300, n_samples = 60, n_modules = 3)
#' net <- rmt_network(mat, verbose = FALSE)
#' plot(net)
#' plot(net, which = "spectrum")
#'
#' @export
plot.RMTnetwork <- function(x, which = "all", ...) {

  which <- match.arg(which,
                     c("all", "spectrum", "mp_fit", "heatmap", "modules"))

  if (which == "all") {
    old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    on.exit(par(old_par))
    .plot_spectrum(x)
    .plot_mp_fit(x)
    .plot_heatmap(x)
    .plot_modules(x)
  } else {
    switch(which,
      spectrum = .plot_spectrum(x),
      mp_fit   = .plot_mp_fit(x),
      heatmap  = .plot_heatmap(x),
      modules  = .plot_modules(x)
    )
  }

  invisible(x)
}

# ── Panel helpers ─────────────────────────────────────────────────────────────

#' @noRd
.plot_spectrum <- function(x) {
  eigs <- x$mp$eigenvalues
  N    <- x$mp$N
  lp   <- x$lambda_plus
  ns   <- x$n_signal

  plot(seq_len(N), eigs,
       type = "l", log = "y",
       xlab = "Eigenvalue rank",
       ylab = "Eigenvalue (log scale)",
       main = "Eigenvalue Spectrum",
       col  = "steelblue", lwd = 1.5)

  abline(h   = lp, col = "firebrick", lty = 2, lwd = 1.5)
  abline(v   = ns, col = "forestgreen", lty = 3, lwd = 1.5)

  legend("topright",
         legend = c("Eigenvalues",
                    sprintf("MP ceiling %.3f", lp),
                    sprintf("Signal cutoff k=%d", ns)),
         col    = c("steelblue", "firebrick", "forestgreen"),
         lty    = c(1, 2, 3), lwd = 1.5, bty = "n", cex = 0.85)
}

#' @noRd
.plot_mp_fit <- function(x) {
  mp    <- x$mp
  sigma2 <- mp$sigma2
  Q      <- mp$Q
  lp     <- mp$lambda_plus
  lm     <- mp$lambda_minus
  eigs   <- mp$eigenvalues

  bulk  <- eigs[eigs <= lp * 1.05]
  grid  <- seq(lm * 0.9, lp * 1.05, length.out = 400)
  pdf_v <- .mp_pdf(grid, sigma2, Q)

  hist(bulk, breaks = 25, freq = FALSE,
       xlab = "Eigenvalue",
       main = "Marchenko-Pastur Fit (bulk)",
       col  = adjustcolor("steelblue", alpha.f = 0.4),
       border = "white")
  lines(grid, pdf_v, col = "firebrick", lwd = 2)
  abline(v = lp, col = "orange", lty = 2, lwd = 1.5)

  legend("topright",
         legend = c("Empirical",
                    sprintf("MP fit (sigma^2=%.3f)", sigma2),
                    sprintf("lambda+ = %.3f", lp)),
         col    = c(adjustcolor("steelblue", 0.6), "firebrick", "orange"),
         lty    = c(NA, 1, 2), lwd = c(NA, 2, 1.5),
         pch    = c(15, NA, NA), bty = "n", cex = 0.85)
}

#' @noRd
.plot_heatmap <- function(x) {
  raw   <- x$filter$raw_corr
  clean <- x$filter$cleaned_corr

  # Simple heatmap of cleaned vs raw difference
  diff_mat <- clean - raw

  # Use image() for a lightweight heatmap
  image(seq_len(nrow(clean)), seq_len(ncol(clean)),
        diff_mat,
        col  = colorRampPalette(c("dodgerblue", "white", "firebrick"))(64),
        xlab = "Gene", ylab = "Gene",
        main = "Noise Removed  (cleaned - raw)")
  box()
}

#' @noRd
.plot_modules <- function(x) {
  mt <- x$module_table
  if (nrow(mt) == 0L) {
    plot.new()
    title(main = "No modules detected")
    return(invisible(NULL))
  }
  cols <- grDevices::hcl.colors(nrow(mt), palette = "Set2")
  barplot(mt$size,
          names.arg = paste0("M", mt$module),
          col       = cols,
          xlab      = "Module",
          ylab      = "Number of genes",
          main      = sprintf("Gene Modules  (n=%d)", nrow(mt)),
          las       = 2)
  abline(h = 0)
}
