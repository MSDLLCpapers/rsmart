# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))


test_that("pstep returns a list with pk, ps, and pfit", {
  p1 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  set.seed(42)
  n <- 100
  a1 <- rbinom(n, size = 1, prob = 0.5)
  df <- data.frame(a1 = a1, x1 = rnorm(n))

  result <- pstep(pmodel = p1, data = df, response = df["a1"], k = 1)

  expect_type(result, "list")
  expect_named(result, c("pk", "ps", "pfit"))
})

test_that("pstep pk values are between 0 and 1", {
  p1 <- modelObj::buildModelObj(
    model = ~ x1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  set.seed(42)
  n <- 200
  x1 <- rnorm(n)
  a1 <- rbinom(n, size = 1, prob = plogis(0.5 * x1))
  df <- data.frame(a1 = a1, x1 = x1)

  result <- pstep(pmodel = p1, data = df, response = df["a1"], k = 1)

  expect_true(all(result$pk >= 0 & result$pk <= 1))
})

test_that("pstep ps values are between 0 and 1", {
  p1 <- modelObj::buildModelObj(
    model = ~ x1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  set.seed(42)
  n <- 200
  x1 <- rnorm(n)
  a1 <- rbinom(n, size = 1, prob = plogis(0.5 * x1))
  df <- data.frame(a1 = a1, x1 = x1)

  result <- pstep(pmodel = p1, data = df, response = df["a1"], k = 1)

  expect_true(all(result$ps >= 0 & result$ps <= 1))
})

test_that("pstep ps equals pk when treatment is 1, and 1-pk when treatment is 0", {
  p1 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  set.seed(42)
  n <- 100
  a1 <- rbinom(n, size = 1, prob = 0.5)
  df <- data.frame(a1 = a1)

  result <- pstep(pmodel = p1, data = df, response = df["a1"], k = 1)

  treated <- df$a1 == 1
  expect_equal(result$ps[treated], result$pk[treated])
  expect_equal(result$ps[!treated], 1 - result$pk[!treated])
})

test_that("pstep works with intercept-only model and equal randomization", {
  p1 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  n <- 1000
  a1 <- c(rep(0, 500), rep(1, 500))
  df <- data.frame(a1 = a1)

  result <- pstep(pmodel = p1, data = df, response = df["a1"], k = 1)

  # With equal randomization and intercept-only model, pk should be 0.5
  expect_equal(mean(result$pk), 0.5, tolerance = 1e-8)
  # All pk values should be identical (intercept-only)
  expect_true(length(unique(result$pk)) == 1)
})

test_that("pstep pk and ps have same length as input data", {
  p1 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  set.seed(42)
  n <- 50
  a1 <- rbinom(n, size = 1, prob = 0.5)
  df <- data.frame(a1 = a1)

  result <- pstep(pmodel = p1, data = df, response = df["a1"], k = 1)

  expect_length(result$pk, n)
  expect_length(result$ps, n)
})

test_that("pstep works for stage k=2", {
  p2 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  set.seed(42)
  n <- 100
  a2 <- rbinom(n, size = 1, prob = 0.5)
  df <- data.frame(a1 = rbinom(n, 1, 0.5), a2 = a2)

  result <- pstep(pmodel = p2, data = df, response = df["a2"], k = 2)

  expect_type(result, "list")
  expect_named(result, c("pk", "ps", "pfit"))
  expect_length(result$pk, n)
  expect_true(all(result$ps >= 0 & result$ps <= 1))
})
