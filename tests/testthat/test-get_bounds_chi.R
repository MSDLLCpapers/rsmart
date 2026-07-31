# ============================================================
# Tests for get_bounds_chi
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# --- Single analysis (K=1) ---

test_that("get_bounds_chi with single analysis and 2 regimes returns correct alpha", {
  result <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                         spend_fn = "OF", corr = diag(2),
                         B = 100001, seed = 42)

  # Type I error should be close to 0.05
  expect_equal(result$typeI, 0.05, tolerance = 0.005)
})

test_that("get_bounds_chi with single analysis and 2 regimes returns bound near chi-squared quantile", {
  result <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                         spend_fn = "OF", corr = diag(2),
                         B = 100001, seed = 42)

  # With 2 independent regimes and contrast matrix of rank 1,
  # the chi-squared statistic has 1 df under the null.
  # The 0.95 quantile of chi-sq(1) is ~3.841
  # use the conversion of the bound to the quantile
  chi_quant = pchisq(result$bound, df = 1, lower.tail = FALSE)
  expect_equal(0.05, chi_quant, tolerance = 0.015)
})

test_that("get_bounds_chi with single analysis returns one bound", {
  result <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                         spend_fn = "OF", corr = diag(2),
                         B = 100001, seed = 42)

  expect_length(result$bound, 1)
})

# --- Single analysis with more regimes ---

test_that("get_bounds_chi with 4 independent regimes has df=3", {
  result <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                         spend_fn = "OF", corr = diag(4),
                         B = 100001, seed = 42)

  # Contrast matrix is 3x4, so chi-sq has 3 df.
  # The 0.95 quantile of chi-sq(3) is ~7.815
  chi_quant = pchisq(result$bound, df = 3, lower.tail = FALSE)
  expect_equal(0.05, chi_quant, tolerance = 0.015)
})

test_that("get_bounds_chi with 4 independent regimes controls alpha", {
  result <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                         spend_fn = "OF", corr = diag(4),
                         B = 100001, seed = 42)

  expect_equal(result$typeI, 0.05, tolerance = 0.005)
})

# --- Two analyses (K=2) ---

test_that("get_bounds_chi with two analyses returns two bounds", {
  result <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                         spend_fn = "OF", corr = diag(2),
                         B = 100001, seed = 42)

  expect_length(result$bound, 2)
})

test_that("get_bounds_chi with two analyses controls overall alpha", {
  result <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                         spend_fn = "OF", corr = diag(2),
                         B = 100001, seed = 42)

  # Overall type I should be close to 0.05
  expect_equal(result$typeI, 0.05, tolerance = 0.005)
})

test_that("get_bounds_chi with OF spending has decreasing bounds", {
  result <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                         spend_fn = "OF", corr = diag(2),
                         B = 100001, seed = 42)

  # O'Brien-Fleming: interim boundary should be more conservative (larger)
  expect_true(result$bound[1] > result$bound[2])
})

test_that("get_bounds_chi with OF spending spends less alpha early", {
  result <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                         spend_fn = "OF", corr = diag(2),
                         B = 100001, seed = 42)

  # OF spends very little alpha at the interim
  expect_true(result$spending[1] < result$spending[2])
})

# --- Pocock spending ---

test_that("get_bounds_chi with Pocock spending controls overall alpha", {
  result <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                         spend_fn = "Pocock", corr = diag(2),
                         B = 100001, seed = 42)

  expect_equal(0.05, result$typeI, tolerance = 0.015)
})

test_that("get_bounds_chi Pocock bounds are more similar than OF bounds", {
  of_result <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                            spend_fn = "OF", corr = diag(2),
                            B = 100001, seed = 42)
  pk_result <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                            spend_fn = "Pocock", corr = diag(2),
                            B = 100001, seed = 42)

  # Pocock bounds should be closer together than OF bounds
  of_diff <- abs(of_result$bound[1] - of_result$bound[2])
  pk_diff <- abs(pk_result$bound[1] - pk_result$bound[2])
  expect_true(pk_diff < of_diff)
})

# --- Reproducibility ---

test_that("get_bounds_chi is reproducible with the same seed", {
  r1 <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                     spend_fn = "OF", corr = diag(2),
                     B = 10001, seed = 123)
  r2 <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                     spend_fn = "OF", corr = diag(2),
                     B = 10001, seed = 123)

  expect_identical(r1$bound, r2$bound)
  expect_identical(r1$typeI, r2$typeI)
})

# --- Different alpha levels ---

test_that("get_bounds_chi with alpha=0.10 has smaller bound than alpha=0.05", {
  r_05 <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                       spend_fn = "OF", corr = diag(2),
                       B = 100001, seed = 42)
  r_10 <- get_bounds_chi(alpha = 0.10, inf_frac = 1,
                       spend_fn = "OF", corr = diag(2),
                       B = 100001, seed = 42)

  # Larger alpha -> smaller (less conservative) boundary
  expect_true(r_10$bound < r_05$bound)
})
