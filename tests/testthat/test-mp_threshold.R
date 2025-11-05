test_that("mp_threshold returns correct structure", {
  set.seed(1)
  mat <- matrix(rnorm(300 * 60), nrow = 300, ncol = 60)
  rownames(mat) <- paste0("g", seq_len(300))

  res <- mp_threshold(mat)

  expect_named(res, c("lambda_plus", "lambda_minus", "sigma2", "Q",
                       "N", "T", "eigenvalues", "eigenvectors",
                       "corr_matrix", "n_signal", "n_noise"))
  expect_equal(res$N, 300L)
  expect_equal(res$T, 60L)
  expect_equal(res$Q, 60 / 300)
  expect_true(res$lambda_plus > res$lambda_minus)
  expect_true(res$sigma2 > 0)
  expect_length(res$eigenvalues, 300L)
  expect_true(all(res$eigenvalues >= 0))
  expect_equal(res$n_signal + res$n_noise, 300L)
})

test_that("pure noise matrix has near-zero signal PCs", {
  set.seed(99)
  # Pure noise — MP law predicts nearly all eigenvalues in bulk
  mat <- matrix(rnorm(500 * 100), nrow = 500, ncol = 100)
  res <- mp_threshold(mat)
  # Allow up to 5% false positives
  expect_lt(res$n_signal / res$N, 0.05)
})

test_that("matrix with embedded modules has detectable signal PCs", {
  mat <- simulate_expression(n_genes = 300, n_samples = 80,
                              n_modules = 3, module_size = 50,
                              signal_strength = 0.8, seed = 7)
  res <- mp_threshold(mat)
  # Should detect at least 3 signal PCs
  expect_gte(res$n_signal, 3L)
})

test_that("mp_threshold rejects bad inputs", {
  expect_error(mp_threshold(matrix(1:9, 3, 3)), "at least 10")
  expect_error(mp_threshold(matrix(c(NA, rnorm(299 * 60 - 1)), 300, 60)),
               "NA")
})

test_that("trace of correlation matrix equals N", {
  set.seed(5)
  mat <- matrix(rnorm(200 * 50), nrow = 200, ncol = 50)
  res <- mp_threshold(mat)
  expect_equal(sum(diag(res$corr_matrix)), res$N, tolerance = 1e-6)
})
