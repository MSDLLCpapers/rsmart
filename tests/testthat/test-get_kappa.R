# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))


# Helper to create test data frames with stage arrival times
make_kappa_df <- function() {
  # 4 individuals in a K=2 trial (columns t1, t2, t3)
  # t1 = enrollment, t2 = stage 2 arrival, t3 = outcome observed
  data.frame(
    t1 = c(1, 5, 10, 20),
    t2 = c(3, 8, 15, 25),
    t3 = c(6, 12, 22, 30)
  )
}

test_that("The length of the returned vector is nrow(df)", {
  df <- make_kappa_df()
  K <- 2
  t_s <- 10

  result <- get_kappa(df, t_s, K)
  expect_length(result, nrow(df))
})

test_that("The values of the returned vector range from 0 (not yet enrolled) to K+1 (outcome observed)", {
  df <- make_kappa_df()
  K <- 2

  # At t_s = 0, no one has enrolled yet — all kappa = 0
  result <- get_kappa(df, t_s = 0, K)
  expect_true(all(result == 0))

  # At t_s = 30, everyone has completed — all kappa = K+1 = 3

  result <- get_kappa(df, t_s = 30, K)
  expect_true(all(result == K + 1))

  # At t_s = 10, individual values should be between 0 and K+1
  result <- get_kappa(df, t_s = 10, K)
  expect_true(all(result >= 0 & result <= K + 1))
})

test_that("get_kappa computes correct stage counts at a given analysis time", {
  df <- make_kappa_df()
  K <- 2

  # At t_s = 10:
  # Individual 1: t1=1<=10, t2=3<=10, t3=6<=10  → kappa = 3
  # Individual 2: t1=5<=10, t2=8<=10, t3=12>10  → kappa = 2
  # Individual 3: t1=10<=10, t2=15>10, t3=22>10 → kappa = 1
  # Individual 4: t1=20>10, t2=25>10, t3=30>10  → kappa = 0
  result <- get_kappa(df, t_s = 10, K)
  expect_equal(result, c(3, 2, 1, 0))
})

test_that("get_kappa throws an error when df is not a data frame", {
  K <- 2
  t_s <- 10

  expect_error(get_kappa(matrix(1:6, nrow = 2), t_s, K), regexp = "data frame")
  expect_error(get_kappa(list(t1 = 1, t2 = 2, t3 = 3), t_s, K), regexp = "data frame")
})

test_that("get_kappa throws an error when t_s is not a non-negative numeric value", {
  df <- make_kappa_df()
  K <- 2

  expect_error(get_kappa(df, t_s = -1, K), regexp = "non-negative numeric")
  expect_error(get_kappa(df, t_s = "10", K), regexp = "non-negative numeric")
  expect_error(get_kappa(df, t_s = c(5, 10), K), regexp = "non-negative numeric")
})

test_that("get_kappa throws an error when K is not a positive integer", {
  df <- make_kappa_df()
  t_s <- 10

  expect_error(get_kappa(df, t_s, K = 0), regexp = "positive integer")
  expect_error(get_kappa(df, t_s, K = -1), regexp = "positive integer")
  expect_error(get_kappa(df, t_s, K = 1.5), regexp = "positive integer")
  expect_error(get_kappa(df, t_s, K = "2"), regexp = "positive integer")
})

