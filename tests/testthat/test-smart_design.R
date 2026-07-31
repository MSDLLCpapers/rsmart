# ============================================================
# Tests for smart_design
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# --- Fixed control, single analysis ---

test_that("smart_design fixed control returns correct structure", {
  result <- smart_design(
    comp_type = "fixed.control",
    test_type = "one-sided",
    alpha = 0.05,
    beta = 0.2,
    delta = c(0, 3),
    n_init = 50,
    inf_frac = 1,
    spend_fn = "OF",
    corr = diag(2),
    variances = c(100, 100)
  )

  expect_type(result, "list")
  expect_true("boundaries" %in% names(result))
  expect_true("sample.size" %in% names(result))
  expect_true("inputs" %in% names(result))
})

test_that("smart_design fixed control achieves desired power", {
  result <- smart_design(
    comp_type = "fixed.control",
    test_type = "one-sided",
    alpha = 0.05,
    beta = 0.2,
    delta = c(0, 3),
    n_init = 50,
    inf_frac = 1,
    spend_fn = "OF",
    corr = diag(2),
    variances = c(100, 100)
  )

  expect_true(result$sample.size$power >= 0.8)
})

test_that("smart_design fixed control echoes inputs", {
  result <- smart_design(
    comp_type = "fixed.control",
    test_type = "one-sided",
    alpha = 0.05,
    beta = 0.2,
    delta = c(0, 3),
    n_init = 50,
    inf_frac = 1,
    spend_fn = "OF",
    corr = diag(2),
    variances = c(100, 100)
  )

  expect_equal(result$inputs$alpha, 0.05)
  expect_equal(result$inputs$beta, 0.2)
  expect_equal(result$inputs$delta, c(0, 3))
  expect_equal(result$inputs$comp_type, "fixed.control")
  expect_equal(result$inputs$test_type, "one-sided")
  expect_equal(result$inputs$spend_fn, "OF")
})

# --- Fixed control, no delta (boundaries only) ---

test_that("smart_design fixed control without delta returns NULL sample.size", {
  result <- smart_design(
    comp_type = "fixed.control",
    test_type = "one-sided",
    alpha = 0.05,
    inf_frac = 1,
    spend_fn = "OF",
    corr = diag(2)
  )

  expect_null(result$sample.size)
  expect_type(result$boundaries, "list")
})

# --- Fixed control, sequential analysis ---

test_that("smart_design fixed control sequential returns two bounds", {
  result <- smart_design(
    comp_type = "fixed.control",
    test_type = "one-sided",
    alpha = 0.05,
    beta = 0.2,
    delta = c(0, 3),
    n_init = 50,
    inf_frac = c(0.5, 1),
    spend_fn = "OF",
    corr = diag(2),
    variances = c(100, 100)
  )

  expect_length(result$boundaries$bound, 2)
  expect_length(result$sample.size$N, 2)
  expect_true(result$sample.size$power >= 0.8)
})

# --- Global chi-squared, single analysis ---

test_that("smart_design global chi-squared returns correct structure", {
  result <- smart_design(
    comp_type = "global.chi.sq",
    alpha = 0.05,
    beta = 0.2,
    delta = c(3, 0),
    n_init = 50,
    inf_frac = 1,
    spend_fn = "OF",
    corr = diag(2),
    variances = c(100, 100),
    Bb = 100001,
    Bc = 10001
  )

  expect_type(result, "list")
  expect_true("boundaries" %in% names(result))
  expect_true("sample.size" %in% names(result))
  expect_true(result$sample.size$power >= 0.8)
})

test_that("smart_design global chi-squared without delta returns NULL sample.size", {
  result <- smart_design(
    comp_type = "global.chi.sq",
    alpha = 0.05,
    inf_frac = 1,
    spend_fn = "OF",
    corr = diag(2),
    Bb = 100001
  )

  expect_null(result$sample.size)
  expect_type(result$boundaries, "list")
})
