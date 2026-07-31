# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))


# Helper to create test data frames
make_test_df <- function(kappa = c(1, 2, 3, 3, 2, 1)) {
  data.frame(kappa = kappa)
}

test_that("get_nu returns a list with nu, ns, and nd", {
  df <- make_test_df()
  K <- 2
  result <- get_nu(df, K)

  expect_type(result, "list")
  expect_named(result, c("nu", "ns", "nd"))
})

test_that("nu has length K+1", {
  df <- make_test_df()
  K <- 2
  result <- get_nu(df, K)

  expect_length(result$nu, K + 1)
})

test_that("ns equals the number of individuals with kappa > 0", {
  # Default data: kappa = c(1, 2, 3, 3, 2, 1) — all enrolled, so ns = 6
  df <- make_test_df()
  K <- 2
  result <- get_nu(df, K)
  expect_equal(result$ns, 6)

  # Mix of enrolled and non-enrolled: kappa = c(0, 1, 2, 3) — ns = 3
  df <- make_test_df(kappa = c(0, 1, 2, 3))
  result <- get_nu(df, K)
  expect_equal(result$ns, 3)
})

test_that("nu values are the proportion of enrolled individuals with kappa >= k", {
  df <- make_test_df()
  K <- 2
  result <- get_nu(df, K)

  # nu[[1]] = 6/6 = 1, nu[[2]] = 4/6, nu[[3]] = 2/6
  expect_equal(result$nu[[1]], 1)
  expect_equal(result$nu[[2]], 4 / 6)
  expect_equal(result$nu[[3]], 2 / 6)
})

test_that("get_nu throws an error when df has no column named 'kappa'", {
  df <- data.frame(other = c(1, 2, 3))
  K <- 2
  expect_error(get_nu(df, K), regexp = "kappa")
})

test_that("get_nu throws an error when ns is not a positive integer", {
  # All kappa values are 0 means no one enrolled (ns = 0), causing division by zero
  df <- make_test_df(kappa = c(0, 0, 0))
  K <- 2

  expect_error(get_nu(df, K), regexp = "ns must be a positive integer")
})

test_that("the nd component of get_nu returns the correct value", {
  df <- make_test_df()
  K <- 2
  result <- get_nu(df, K)

  # nd = sum(kappa >= K+1) / sum(kappa >= K) = sum(kappa >= 3) / sum(kappa >= 2) = 2/4
  expect_equal(result$nd, 0.5)
})

test_that("get_nu throws an error when there would be division by zero in nd calculation", {
  # Everyone is enrolled (kappa > 0) so ns check passes, but no one reaches
  # stage K, so sum(kappa >= K) = 0 and nd = 0/0
  df <- make_test_df(kappa = c(1, 1, 1))
  K <- 3

  expect_error(get_nu(df, K), regexp = "nd denominator")
})

test_that("get_nu handles the boundary cases of everyone finishes & no one finishes", {

  # case 1: everyone completes the full trial (kappa = K + 1 for all)
  df <- make_test_df(kappa = c(3, 3, 3))
  K <- 2
  result <- get_nu(df, K)

  expect_equal(result$nu[[1]], 1)
  expect_equal(result$nu[[2]], 1)
  expect_equal(result$nu[[3]], 1)

  # case 2: everyone enrolled but no one progresses past stage 1
  # nu[[1]] = 4/4 = 1, nu[[2]] = 2/4, nu[[3]] = 0/4 = 0
  df <- make_test_df(kappa = c(1, 2, 2, 1))
  K <- 2
  result <- get_nu(df, K)

  expect_equal(result$nu[[1]], 1)
  expect_equal(result$nu[[2]], 2 / 4)
  expect_equal(result$nu[[3]], 0)
})

test_that("get_nu works correctly for a single-stage trial (K = 1)", {
  # kappa = 1 means enrolled but not finished; kappa = 2 means outcome observed
  df <- make_test_df(kappa = c(1, 2, 2, 1, 2))
  K <- 1
  result <- get_nu(df, K)

  # nu should have K+1 = 2 elements
  expect_length(result$nu, 2)

  # nu[[1]] = 5/5 = 1 (all enrolled)
  # nu[[2]] = 3/5 (three have kappa >= 2, i.e. outcome observed)
  expect_equal(result$nu[[1]], 1)
  expect_equal(result$nu[[2]], 3 / 5)

  # ns = 5

  expect_equal(result$ns, 5)

  # nd = sum(kappa >= 2) / sum(kappa >= 1) = 3/5
  expect_equal(result$nd, 3 / 5)
})
