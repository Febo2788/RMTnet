test_that("rmt_network returns RMTnetwork object", {
  mat <- simulate_expression(n_genes = 200, n_samples = 50,
                              n_modules = 2, seed = 1)
  net <- rmt_network(mat, min_module_size = 10, verbose = FALSE)

  expect_s3_class(net, "RMTnetwork")
  expect_true(is.numeric(net$entropy))
  expect_true(net$entropy > 0)
  expect_true(net$entropy <= net$entropy_max)
  expect_equal(nrow(net$cleaned_corr), 200L)
  expect_equal(ncol(net$cleaned_corr), 200L)
  expect_true(all(diag(net$cleaned_corr) == 1))
})

test_that("rmt_network detects embedded modules", {
  mat <- simulate_expression(n_genes = 300, n_samples = 80,
                              n_modules = 3, module_size = 50,
                              signal_strength = 0.8, seed = 42)
  net <- rmt_network(mat, min_module_size = 20, verbose = FALSE)

  # Should find at least 2 of the 3 embedded modules
  expect_gte(nrow(net$module_table), 2L)
})

test_that("print.RMTnetwork works without error", {
  mat <- simulate_expression(n_genes = 150, n_samples = 40,
                              n_modules = 2, seed = 3)
  net <- rmt_network(mat, min_module_size = 10, verbose = FALSE)
  expect_output(print(net), "RMTnet")
})
