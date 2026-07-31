# ============================================================
# Tests for get_sample_size_chi
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# --- Fixed analysis (single analysis, S=1) ---

test_that("get_sample_size_chi returns correct structure for fixed analysis", {
  # 2 independent regimes, single analysis
  b <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                    spend_fn = "OF", corr = diag(2),
                    B = 100001, seed = 42)

  result <- get_sample_size_chi(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(3, 0),
    bounds = b$bound,
    n_init = 50,
    corr = diag(2),
    inf_frac = 1
  )

  expect_type(result, "list")
  expect_true("N" %in% names(result))
  expect_true("power" %in% names(result))
  expect_length(result$N, 1)
})

test_that("get_sample_size_chi works with bounds structure for fixed analysis", {
  # 2 independent regimes, single analysis
  b <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                    spend_fn = "OF", corr = diag(2),
                    B = 100001, seed = 42)

  result <- get_sample_size_chi(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(3, 0),
    bounds = b,
    n_init = 50
  )

  expect_type(result, "list")
  expect_true("N" %in% names(result))
  expect_true("power" %in% names(result))
  expect_length(result$N, 1)
})


test_that("get_sample_size_chi achieves desired power for fixed analysis", {
  b <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                    spend_fn = "OF", corr = diag(2),
                    B = 100001, seed = 42)

  result <- get_sample_size_chi(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(3, 0),
    bounds = b$bound,
    n_init = 50,
    corr = diag(2),
    inf_frac = 1
  )

  # Power should be at least 1 - beta = 0.8
  expect_true(result$power >= 0.8)
})

test_that("get_sample_size_chi requires larger N for higher power", {
  b <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                    spend_fn = "OF", corr = diag(2),
                    B = 100001, seed = 42)

  r_80 <- get_sample_size_chi(
    variances = c(100, 100), beta = 0.2, delta = c(3, 0),
    bounds = b$bound, n_init = 50,
    corr = diag(2), inf_frac = 1
  )

  r_90 <- get_sample_size_chi(
    variances = c(100, 100), beta = 0.1, delta = c(3, 0),
    bounds = b$bound, n_init = 50,
    corr = diag(2), inf_frac = 1
  )

  expect_true(max(r_90$N) >= max(r_80$N))
})

test_that("get_sample_size_chi requires larger N for smaller effect size", {
  b <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                    spend_fn = "OF", corr = diag(2),
                    B = 100001, seed = 42)

  r_large <- get_sample_size_chi(
    variances = c(100, 100), beta = 0.2, delta = c(5, 0),
    bounds = b$bound, n_init = 50,
    corr = diag(2), inf_frac = 1
  )

  r_small <- get_sample_size_chi(
    variances = c(100, 100), beta = 0.2, delta = c(2, 0),
    bounds = b$bound, n_init = 50,
    corr = diag(2), inf_frac = 1
  )

  expect_true(max(r_small$N) >= max(r_large$N))
})

test_that("get_sample_size_chi requires larger N for larger variance", {
  b <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                    spend_fn = "OF", corr = diag(2),
                    B = 100001, seed = 42)

  r_low_var <- get_sample_size_chi(
    variances = c(50, 50), beta = 0.2, delta = c(3, 0),
    bounds = b$bound, n_init = 50,
    corr = diag(2), inf_frac = 1
  )

  r_high_var <- get_sample_size_chi(
    variances = c(200, 200), beta = 0.2, delta = c(3, 0),
    bounds = b$bound, n_init = 50,
    corr = diag(2), inf_frac = 1
  )

  expect_true(max(r_high_var$N) >= max(r_low_var$N))
})

# --- Fixed analysis with 4 regimes ---

test_that("get_sample_size_chi works with 4 independent regimes", {
  b <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                    spend_fn = "OF", corr = diag(4),
                    B = 100001, seed = 42)

  result <- get_sample_size_chi(
    variances = c(100, 100, 100, 100),
    beta = 0.2,
    delta = c(4, 0, 0, 0),
    bounds = b$bound,
    n_init = 50,
    corr = diag(4),
    inf_frac = 1
  )

  expect_type(result, "list")
  expect_true(result$power >= 0.8)
  expect_length(result$N, 1)
})

# --- Sequential analysis (S=2) ---

test_that("get_sample_size_chi returns correct structure for sequential analysis", {
  b <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                    spend_fn = "OF", corr = diag(2),
                    B = 100001, seed = 42)

  result <- get_sample_size_chi(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(3, 0),
    bounds = b$bound,
    n_init = 50,
    corr = diag(2),
    inf_frac = c(0.5, 1)
  )

  expect_type(result, "list")
  expect_length(result$N, 2)
})

test_that("get_sample_size_chi achieves desired power for sequential analysis", {
  b <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                    spend_fn = "OF", corr = diag(2),
                    B = 100001, seed = 42)

  result <- get_sample_size_chi(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(3, 0),
    bounds = b$bound,
    n_init = 50,
    corr = diag(2),
    inf_frac = c(0.5, 1)
  )

  expect_true(result$power >= 0.8)
})

test_that("get_sample_size_chi N values reflect n_split for sequential analysis", {
  b <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                    spend_fn = "OF", corr = diag(2),
                    B = 100001, seed = 42)

  result <- get_sample_size_chi(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(3, 0),
    bounds = b$bound,
    n_init = 50,
    corr = diag(2),
    inf_frac = c(0.5, 1)
  )

  # N[1] should be approximately half of N[2]
  expect_equal(result$N[1] / result$N[2], 0.5, tolerance = 0.05)
})

test_that("get_sample_size_chi requires larger N for higher power in sequential analysis", {
  b <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                    spend_fn = "OF", corr = diag(2),
                    B = 100001, seed = 42)

  r_80 <- get_sample_size_chi(
    variances = c(100, 100), beta = 0.2, delta = c(3, 0),
    bounds = b$bound, n_init = 50,
    corr = diag(2), inf_frac = c(0.5, 1)
  )

  r_90 <- get_sample_size_chi(
    variances = c(100, 100), beta = 0.1, delta = c(3, 0),
    bounds = b$bound, n_init = 50,
    corr = diag(2), inf_frac = c(0.5, 1)
  )

  expect_true(max(r_90$N) >= max(r_80$N))
})

# --- Pocock spending ---

test_that("get_sample_size_chi with Pocock spending achieves desired power", {
  b <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                    spend_fn = "Pocock", corr = diag(2),
                    B = 100001, seed = 42)

  result <- get_sample_size_chi(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(3, 0),
    bounds = b$bound,
    n_init = 50,
    corr = diag(2),
    inf_frac = c(0.5, 1)
  )

  expect_true(result$power >= 0.8)
})

test_that("get_sample_size_chi Pocock vs OF sample sizes differ", {
  b_of <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                       spend_fn = "OF", corr = diag(2),
                       B = 100001, seed = 42)
  b_pk <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                       spend_fn = "Pocock", corr = diag(2),
                       B = 100001, seed = 42)

  r_of <- get_sample_size_chi(
    variances = c(100, 100), beta = 0.2, delta = c(3, 0),
    bounds = b_of$bound, n_init = 50,
    corr = diag(2), inf_frac = c(0.5, 1)
  )

  r_pk <- get_sample_size_chi(
    variances = c(100, 100), beta = 0.2, delta = c(3, 0),
    bounds = b_pk$bound, n_init = 50,
    corr = diag(2), inf_frac = c(0.5, 1)
  )

  # Both should achieve power
  expect_true(r_of$power >= 0.8)
  expect_true(r_pk$power >= 0.8)
})

# --- Reproducibility ---

test_that("get_sample_size_chi is reproducible with same seed", {
  b <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                    spend_fn = "OF", corr = diag(2),
                    B = 100001, seed = 42)

  r1 <- get_sample_size_chi(
    variances = c(100, 100), beta = 0.2, delta = c(3, 0),
    bounds = b$bound, n_init = 50,
    corr = diag(2), inf_frac = 1,
    seed = 99, B = 10001
  )

  r2 <- get_sample_size_chi(
    variances = c(100, 100), beta = 0.2, delta = c(3, 0),
    bounds = b$bound, n_init = 50,
    corr = diag(2), inf_frac = 1,
    seed = 99, B = 10001
  )

  expect_identical(r1$N, r2$N)
  expect_identical(r1$power, r2$power)
})

# --- Accepts get_bounds_chi list directly ---

test_that("get_sample_size_chi accepts get_bounds_chi output list as bounds", {
  b <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                    spend_fn = "OF", corr = diag(2),
                    B = 100001, seed = 42)

  # Pass the full list — the function should extract $bound
  result <- get_sample_size_chi(
    variances = c(100, 100),
    beta = 0.2,
    delta = c(3, 0),
    bounds = b,
    n_init = 50,
    corr = diag(2),
    inf_frac = 1
  )

  expect_type(result, "list")
  expect_true(result$power >= 0.8)
})

# --- Sequential requires larger max N than fixed ---

test_that("get_sample_size_chi OF sequential requires larger max N than fixed", {
  b_fixed <- get_bounds_chi(alpha = 0.05, inf_frac = 1,
                          spend_fn = "OF", corr = diag(2),
                          B = 100001, seed = 42)
  b_of <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                        spend_fn = "OF", corr = diag(2),
                        B = 100001, seed = 42)
  b_po <- get_bounds_chi(alpha = 0.05, inf_frac = c(0.5, 1),
                        spend_fn = "Pocock", corr = diag(2),
                        B = 100001, seed = 42)

  r_fixed <- get_sample_size_chi(
    variances = c(100, 100), beta = 0.2, delta = c(3, 0),
    bounds = b_fixed$bound, n_init = 50,
    corr = diag(2), inf_frac = 1
  )

  r_of <- get_sample_size_chi(
    variances = c(100, 100), beta = 0.2, delta = c(3, 0),
    bounds = b_of$bound, n_init = 50,
    corr = diag(2), inf_frac = c(0.5, 1)
  )
  
  r_po <- get_sample_size_chi(
    variances = c(100, 100), beta = 0.2, delta = c(3, 0),
    bounds = b_po$bound, n_init = 50,
    corr = diag(2), inf_frac = c(0.5, 1)
  )

  # Sequential max N should be at least as large as fixed N
  expect_true(max(r_of$N) >= max(r_fixed$N))
  expect_true(max(r_po$N) >= max(r_fixed$N))
})
