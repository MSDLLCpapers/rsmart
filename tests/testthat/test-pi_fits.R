# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))


test_that("pi_fits returns a list with ps and p_fits", {
  p1 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  set.seed(42)
  n <- 100
  df <- data.frame(
    a1 = rbinom(n, size = 1, prob = 0.5),
    kappa = rep(1, n)
  )

  result <- pi_fits(df = df, p_list = list(p1))

  expect_type(result, "list")
  expect_named(result, c("ps", "p_fits"))
})

test_that("pi_fits ps has correct dimensions and column names", {
  p1 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )
  p2 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  set.seed(42)
  n <- 100
  df <- data.frame(
    a1 = rbinom(n, size = 1, prob = 0.5),
    a2 = rbinom(n, size = 1, prob = 0.5),
    kappa = c(rep(1, 50), rep(2, 50))
  )

  result <- pi_fits(df = df, p_list = list(p1, p2))

  expect_equal(nrow(result$ps), n)
  expect_equal(ncol(result$ps), 2)
  expect_equal(colnames(result$ps), c("pi1", "pi2"))
})

test_that("pi_fits returns one fitted model per stage", {
  p1 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )
  p2 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  set.seed(42)
  n <- 100
  df <- data.frame(
    a1 = rbinom(n, size = 1, prob = 0.5),
    a2 = rbinom(n, size = 1, prob = 0.5),
    kappa = rep(3, n)
  )

  result <- pi_fits(df = df, p_list = list(p1, p2))

  expect_length(result$p_fits, 2)
})

test_that("pi_fits assigns 99 for individuals who have not reached a stage", {
  p1 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )
  p2 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  # first 50 have kappa=1 (only reached stage 1), last 50 have kappa=2
  df <- data.frame(
    a1 = c(rep(0, 50), rep(1, 50), rep(0, 50), rep(1, 50)),
    a2 = c(rep(0, 50), rep(1, 50), rep(0, 50), rep(1, 50)),
    kappa = c(rep(c(1, 2), 100))
  )

  result <- pi_fits(df = df, p_list = list(p1, p2))

  # all individuals should have a pi1 estimate (everyone reached stage 1)
  expect_true(all(result$ps$pi1 != 99))

  # individuals with kappa=1 should have 99 for pi2
  expect_true(all(result$ps$pi2[df$kappa == 1] == 99))
  # individuals with kappa=2 should have a real estimate for pi2
  expect_true(all(result$ps$pi2[df$kappa >= 2] != 99))
  # all other individuals should have pi2 estimates of 0.5
  expect_true(all(result$ps$pi2[df$kappa >= 2] == 0.5))
})

test_that("pi_fits ps values are between 0 and 1 for estimated individuals", {
  p1 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )
  p2 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  set.seed(42)
  n <- 200
  df <- data.frame(
    a1 = rbinom(n, size = 1, prob = 0.5),
    a2 = rbinom(n, size = 1, prob = 0.5),
    kappa = sample(1:3, n, replace = TRUE)
  )

  result <- pi_fits(df = df, p_list = list(p1, p2))

  # check only estimated values (not the 99 placeholders)
  pi1_est <- result$ps$pi1[result$ps$pi1 != 99]
  pi2_est <- result$ps$pi2[result$ps$pi2 != 99]
  expect_true(all(pi1_est >= 0 & pi1_est <= 1))
  expect_true(all(pi2_est >= 0 & pi2_est <= 1))
})

test_that("pi_fits with K=2 and equally divided treatment groups gives ps 0.5", {
  p1 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )
  p2 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  # Create data with exactly equal groups: (0,0), (0,1), (1,0), (1,1)
  # each with 50 individuals, all reaching stage 2
  n_per_group <- 50
  df <- data.frame(
    a1 = rep(c(0, 0, 1, 1), each = n_per_group),
    a2 = rep(c(0, 1, 0, 1), each = n_per_group),
    kappa = rep(3, 4 * n_per_group)
  )

  result <- pi_fits(df = df, p_list = list(p1, p2))

  # With exactly equal treatment splits and intercept-only models,
  # the estimated probability of treatment=1 should be exactly 0.5
  # so ps should be 0.5 for everyone
  expect_equal(unique(result$ps$pi1), 0.5)
  expect_equal(unique(result$ps$pi2), 0.5)

  # No placeholders since everyone reached both stages
  expect_true(all(result$ps$pi1 != 99))
  expect_true(all(result$ps$pi2 != 99))

  # With intercept-only models on balanced data, pk = 0.5 for both stages.
  # ps = pk when a=1, ps = 1-pk when a=0. Since pk=0.5, ps=0.5 for all.
  expect_true(all(abs(result$ps$pi1 - 0.5) < 1e-10))
  expect_true(all(abs(result$ps$pi2 - 0.5) < 1e-10))
})

test_that("pi_fits with K=2 equal groups and unbalanced kappa handles subsetting", {
  p1 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )
  p2 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  # 200 total: 100 reach only stage 1 (kappa=1), 100 reach stage 2 (kappa=2)
  # Within each kappa group, a1 and a2 are equally split
  n_half <- 100
  df <- data.frame(
    a1 = c(rep(0, n_half / 2), rep(1, n_half / 2),
            rep(0, n_half / 2), rep(1, n_half / 2)),
    a2 = c(rep(0, n_half / 4), rep(1, n_half / 4),
            rep(0, n_half / 4), rep(1, n_half / 4),
            rep(0, n_half / 4), rep(1, n_half / 4),
            rep(0, n_half / 4), rep(1, n_half / 4)),
    kappa = c(rep(1, n_half), rep(2, n_half))
  )

  result <- pi_fits(df = df, p_list = list(p1, p2))

  # Stage 1: all 200 individuals used, a1 equally split -> pi1 = 0.5
  expect_true(all(abs(result$ps$pi1 - 0.5) < 1e-10))

  # Stage 2: only kappa>=2 individuals (last 100) used
  # Those with kappa=1 should have placeholder 99

  expect_true(all(result$ps$pi2[df$kappa == 1] == 99))
  # Those with kappa>=2 should have ps = 0.5
  expect_true(all(abs(result$ps$pi2[df$kappa >= 2] - 0.5) < 1e-10))
})

test_that("pi_fits works with a single stage (K=1)", {
  p1 <- modelObj::buildModelObj(
    model = ~ 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  n <- 80
  df <- data.frame(
    a1 = c(rep(0, 40), rep(1, 40)),
    kappa = rep(1, n)
  )

  result <- pi_fits(df = df, p_list = list(p1))

  expect_equal(ncol(result$ps), 1)
  expect_equal(colnames(result$ps), "pi1")
  expect_length(result$p_fits, 1)
  expect_true(all(abs(result$ps$pi1 - 0.5) < 1e-10))
})
