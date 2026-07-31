# ============================================================
# Tests for smart_design and get_bounds
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# Test get bounds with single analysis 

test_that("get_bounds with single analysis and 2 independent regimes returns correct alpha", {
  result <- get_bounds(alpha = 0.05, inf_frac = 1,
                      spend_fn = "OF", corr = diag(2), test_type = "one-sided",
                      lambda = 0.1, tol = 1e-6, max_iter = 1000)

  # The total type I error (spending) should be within 0.00001 of 0.05

  expect_equal(result$spending[1], 0.05, tolerance = 1e-4)
})

test_that("get_bounds with single analysis and 2 independent regimes returns correct bound", {
  result <- get_bounds(alpha = 0.05, inf_frac = 1,
                      spend_fn = "OF", corr = diag(2), test_type = "one-sided",
                      lambda = 0.1, tol = 1e-6, max_iter = 1000)

  # With 2 independent regimes and FWER = 0.05, the bound should be
  # close to the Bonferroni-adjusted z-critical value qnorm(1 - 0.025)
  z_025 <- qnorm(1 - 0.025)
  expect_equal(result$bounds, z_025, tolerance = 0.01)
})

test_that("get_bounds converges with single analysis", {
  result <- get_bounds(alpha = 0.05, inf_frac = 1,
                      spend_fn = "OF", corr = diag(2), test_type = "one-sided",
                      lambda = 0.1, tol = 1e-6, max_iter = 1000)

  expect_true(result$convergence)
})

# Test get bounds with sequential results

test_that("get_bounds with two analyses and 2 independent regimes controls alpha", {
  # corr is 4x4 = 2 regimes x 2 analyses
  result <- get_bounds(alpha = 0.05, inf_frac = c(0.7, 1.0),
                      spend_fn = "OF", corr = diag(2), test_type = "one-sided",
                      lambda = 0.1, tol = 1e-6, max_iter = 1000)

  # Cumulative spending at the final analysis should be within tolerance of 0.05
  expect_equal(result$spending[2], 0.05, tolerance = 1e-4)
})

test_that("get_bounds with two analyses returns two correct bounds", {
  result <- get_bounds(alpha = 0.05, inf_frac = c(0.7, 1.0),
                      spend_fn = "OF", corr = diag(2), test_type = "one-sided",
                      lambda = 0.1, tol = 1e-6, max_iter = 1000)

  expect_length(result$bound, 2)
  expect_true(result$bound[1] > result$bound[2])
  z_1 <- qnorm(1 - result$spending[1]/2)
  expect_equal(result$bound[1], z_1, tolerance = 0.01)
})

test_that("get_bounds with two analyses converges at both stages", {
  result <- get_bounds(alpha = 0.05, inf_frac = c(0.7, 1.0),
                      spend_fn = "OF", corr = diag(2), test_type = "one-sided",
                      lambda = 0.1, tol = 1e-6, max_iter = 1000)

  expect_true(all(result$convergence))
})


# Tests for two-sided analyses
test_that("get_bounds with single analysis and 2 independent regimes returns correct bound", {
  result <- get_bounds(alpha = 0.05, inf_frac = 1,
                      spend_fn = "OF", corr = diag(2), test_type = "two-sided",
                      lambda = 0.1, tol = 1e-6, max_iter = 1000)

  # With 2 independent regimes and FWER = 0.05, the bound should be
  # close to the Bonferroni-adjusted z-critical value qnorm(1 - 0.025/2)
  z_crit <- qnorm(1 - 0.05/4)
  expect_equal(result$bounds, z_crit, tolerance = 0.01)
})

test_that("get_bounds with two analyses and 2 independent regimes controls alpha", {
  # corr is 4x4 = 2 regimes x 2 analyses
  result <- get_bounds(alpha = 0.05, inf_frac = c(0.7, 1.0),
                      spend_fn = "OF", corr = diag(2), test_type = "two-sided",
                      lambda = 0.1, tol = 1e-6, max_iter = 1000)

  # Cumulative spending at the final analysis should be within tolerance of 0.05
  expect_equal(result$spending[2], 0.05, tolerance = 1e-4)
})

test_that("get_bounds with two analyses returns two correct bounds", {
  result <- get_bounds(alpha = 0.05, inf_frac = c(0.7, 1.0),
                      spend_fn = "OF", corr = diag(2), test_type = "two-sided",
                      lambda = 0.1, tol = 1e-6, max_iter = 1000)

  expect_length(result$bound, 2)
  expect_true(result$bound[1] > result$bound[2])
  z_1 <- qnorm(1 - result$spending[1]/4)
  expect_equal(result$bound[1], z_1, tolerance = 0.01)
})

test_that("get_bounds with two analyses converges at both stages", {
  result <- get_bounds(alpha = 0.05, inf_frac = c(0.7, 1.0),
                      spend_fn = "OF", corr = diag(2), test_type = "two-sided",
                      lambda = 0.1, tol = 1e-6, max_iter = 1000)

  expect_true(all(result$convergence))
})
