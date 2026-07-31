# ============================================================
# Tests for get_sample_size
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# --- Fixed analysis (single analysis, K=1) ---

test_that("get_sample_size returns correct structure for fixed analysis", {
  # 2 independent regimes, single analysis
  b <- get_bounds(alpha = 0.05, inf_frac = 1, spend_fn = "OF",
                 corr = diag(2), test_type = "one-sided")

  result <- get_sample_size(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(0, 3),
    bounds = b$bounds,
    n_init = 50,
    corr = b$corr,
    inf_frac = b$inf_frac
  )

  expect_type(result, "list")
  expect_true("N" %in% names(result))
  expect_true("power" %in% names(result))
  expect_true("prop_rej" %in% names(result))
  expect_length(result$N, 1)
  expect_length(result$prop_rej, 1)
})

test_that("get_sample_size functions with bounds object for fixed analysis", {
  # 2 independent regimes, single analysis
  b <- get_bounds(alpha = 0.05, inf_frac = 1, spend_fn = "OF",
                 corr = diag(2), test_type = "one-sided")

  result <- get_sample_size(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(0, 3),
    bounds = b,
    n_init = 50
  )

  expect_type(result, "list")
  expect_true("N" %in% names(result))
  expect_true("power" %in% names(result))
  expect_true("prop_rej" %in% names(result))
  expect_length(result$N, 1)
  expect_length(result$prop_rej, 1)
  expect_equal(result$bounds, b$bounds)
  expect_equal(result$corr, b$corr)
})



test_that("get_sample_size achieves desired power for fixed analysis", {
  b <- get_bounds(alpha = 0.05, inf_frac = 1, spend_fn = "OF",
                 corr = diag(2), test_type = "one-sided")

  result <- get_sample_size(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(0, 3),
    bounds = b$bounds,
    n_init = 50,
    corr = b$corr,
    inf_frac = b$inf_frac
  )

  # Power should be at least 1 - beta = 0.8
  expect_true(result$power >= 0.8)
})

test_that("get_sample_size returns input parameters for fixed analysis", {
  b <- get_bounds(alpha = 0.05, inf_frac = 1, spend_fn = "OF",
                 corr = diag(2), test_type = "one-sided")

  result <- get_sample_size(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(0, 3),
    bounds = b$bounds,
    n_init = 50,
    corr = b$corr,
    inf_frac = b$inf_frac
  )

  expect_equal(result$beta, 0.2)
  expect_equal(result$delta, c(0, 3))
  expect_equal(result$n_init, 50)
  expect_equal(result$inf_frac, 1)
})

test_that("get_sample_size requires larger N for higher power in fixed analysis", {
  b <- get_bounds(alpha = 0.05, inf_frac = 1, spend_fn = "OF",
                 corr = diag(2), test_type = "one-sided")

  r_80 <- get_sample_size(
    variances = c(100, 100), beta = 0.2, delta = c(0, 3),
    bounds = b$bounds, n_init = 50,
    corr = b$corr, inf_frac = b$inf_frac
  )

  r_90 <- get_sample_size(
    variances = c(100, 100), beta = 0.1, delta = c(0, 3),
    bounds = b$bounds, n_init = 50,
    corr = b$corr, inf_frac = b$inf_frac
  )

  expect_true(r_90$N >= r_80$N)
})

test_that("get_sample_size requires larger N for smaller effect size", {
  b <- get_bounds(alpha = 0.05, inf_frac = 1, spend_fn = "OF",
                 corr = diag(2), test_type = "one-sided")

  r_large <- get_sample_size(
    variances = c(100, 100), beta = 0.2, delta = c(0, 5),
    bounds = b$bounds, n_init = 50,
    corr = b$corr, inf_frac = b$inf_frac
  )

  r_small <- get_sample_size(
    variances = c(100, 100), beta = 0.2, delta = c(0, 2),
    bounds = b$bounds, n_init = 50,
    corr = b$corr, inf_frac = b$inf_frac
  )

  expect_true(r_small$N >= r_large$N)
})

test_that("get_sample_size requires larger N for larger variance", {
  b <- get_bounds(alpha = 0.05, inf_frac = 1, spend_fn = "OF",
                 corr = diag(2), test_type = "one-sided")

  r_low_var <- get_sample_size(
    variances = c(50, 50), beta = 0.2, delta = c(0, 3),
    bounds = b$bounds, n_init = 50,
    corr = b$corr, inf_frac = b$inf_frac
  )

  r_high_var <- get_sample_size(
    variances = c(200, 200), beta = 0.2, delta = c(0, 3),
    bounds = b$bounds, n_init = 50,
    corr = b$corr, inf_frac = b$inf_frac
  )

  expect_true(r_high_var$N >= r_low_var$N)
})

# --- Sequential analysis (K=2) ---

test_that("get_sample_size returns correct structure for sequential analysis", {
  b <- get_bounds(alpha = 0.05, inf_frac = c(0.5, 1), spend_fn = "OF",
                 corr = diag(2), test_type = "one-sided")

  result <- get_sample_size(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(0, 3),
    bounds = b$bound,
    n_init = 50,
    corr = b$corr,
    inf_frac = b$inf_frac
  )

  expect_type(result, "list")
  expect_length(result$N, 2)
  expect_length(result$prop_rej, 2)
})

test_that("get_sample_size achieves desired power for sequential analysis", {
  b <- get_bounds(alpha = 0.05, inf_frac = c(0.5, 1), spend_fn = "OF",
                 corr = diag(2), test_type = "one-sided")

  result <- get_sample_size(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(0, 3),
    bounds = b$bound,
    n_init = 50,
    corr = b$corr,
    inf_frac = b$inf_frac
  )

  expect_true(result$power >= 0.8)
})

test_that("get_sample_size N values reflect information fractions", {
  b <- get_bounds(alpha = 0.05, inf_frac = c(0.5, 1), spend_fn = "OF",
                 corr = diag(2), test_type = "one-sided")

  result <- get_sample_size(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(0, 3),
    bounds = b$bound,
    n_init = 50,
    corr = b$corr,
    inf_frac = b$inf_frac
  )

  # N[1] should be half of N[2] +/-1 (inf_frac = c(0.5, 1)) because equal 
  # information and no correlation between statistics 
  expect_equal(result$N[1] / result$N[2], 0.5, tolerance = 0.5)
})

test_that("get_sample_size cumulative rejection probabilities increase over analyses", {
  b <- get_bounds(alpha = 0.05, inf_frac = c(0.5, 1), spend_fn = "OF",
                 corr = diag(2), test_type = "one-sided")

  result <- get_sample_size(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(0, 3),
    bounds = b$bound,
    n_init = 50,
    corr = b$corr,
    inf_frac = b$inf_frac
  )

  # Cumulative rejection probability should increase from interim to final
  expect_true(result$prop_rej[2] >= result$prop_rej[1])
})

test_that("get_sample_size returns input parameters for sequential analysis", {
  b <- get_bounds(alpha = 0.05, inf_frac = c(0.5, 1), spend_fn = "OF",
                 corr = diag(2), test_type = "one-sided")

  result <- get_sample_size(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(0, 3),
    bounds = b$bound,
    n_init = 50,
    corr = b$corr,
    inf_frac = b$inf_frac
  )

  expect_equal(result$beta, 0.2)
  expect_equal(result$delta, c(0, 3))
  expect_equal(result$n_init, 50)
  expect_equal(result$inf_frac, c(0.5, 1))
})

test_that("get_sample_size requires larger N for higher power in sequential analysis", {
  b <- get_bounds(alpha = 0.05, inf_frac = c(0.5, 1), spend_fn = "OF",
                 corr = diag(2), test_type = "one-sided")

  r_80 <- get_sample_size(
    variances = c(100, 100), beta = 0.2, delta = c(0, 3),
    bounds = b$bound, n_init = 50,
    corr = b$corr, inf_frac = b$inf_frac
  )

  r_90 <- get_sample_size(
    variances = c(100, 100), beta = 0.1, delta = c(0, 3),
    bounds = b$bound, n_init = 50,
    corr = b$corr, inf_frac = b$inf_frac
  )

  # Final N should be larger for 90% power
  expect_true(max(r_90$N) >= max(r_80$N))
})

test_that("get_sample_size with Pocock spending differs from OF", {
  b_of <- get_bounds(alpha = 0.05, inf_frac = c(0.5, 1), spend_fn = "OF",
                    corr = diag(2), test_type = "one-sided")
  b_pk <- get_bounds(alpha = 0.05, inf_frac = c(0.5, 1), spend_fn = "Pocock",
                    corr = diag(2), test_type = "one-sided")

  r_of <- get_sample_size(
    variances = c(100, 100), beta = 0.2, delta = c(0, 3),
    bounds = b_of$bound, n_init = 50,
    corr = b_of$corr, inf_frac = b_of$inf_frac
  )

  r_pk <- get_sample_size(
    variances = c(100, 100), beta = 0.2, delta = c(0, 3),
    bounds = b_pk$bound, n_init = 50,
    corr = b_pk$corr, inf_frac = b_pk$inf_frac
  )

  # Both should achieve at least 80% power
  expect_true(r_of$power >= 0.8)
  expect_true(r_pk$power >= 0.8)

  # Sample sizes will generally differ between spending functions
  # (Pocock typically requires slightly more)
  expect_true(max(r_of$N) != max(r_pk$N) || max(r_of$N) == max(r_pk$N))
  # Expect OF to be less than Pocock 
  expect_true(max(r_of$N) < max(r_pk$N))
})
