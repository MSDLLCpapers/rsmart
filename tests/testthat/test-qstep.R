# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))


test_that("qstep returns hats_mod, hats_unmod, and qfit", {
  qmodel <- modelObj::buildModelObj(
    model = ~ x1 + a1 + x1:a1,
    solver.method = "lm",
    predict.method = "predict.lm"
  )

  set.seed(1)
  n <- 120
  dat <- data.frame(
    x1 = rnorm(n),
    a1 = rbinom(n, 1, 0.5)
  )
  dat$y <- with(dat, 0.5 + 0.8 * x1 + 1.2 * a1 + 0.7 * x1 * a1 + rnorm(n, sd = 0.1))

  out <- qstep(
    qmodel = qmodel,
    data = dat,
    response = dat["y"],
    newdata = dat,
    regime = rep(1L, n),
    txName = "a1"
  )

  expect_type(out, "list")
  expect_named(out, c("hats_mod", "hats_unmod", "qfit"))
})

test_that("qstep predictions have expected lengths", {
  qmodel <- modelObj::buildModelObj(
    model = ~ x1 + a1,
    solver.method = "lm",
    predict.method = "predict.lm"
  )

  set.seed(2)
  n_fit <- 100
  n_new <- 40
  fit_dat <- data.frame(
    x1 = rnorm(n_fit),
    a1 = rbinom(n_fit, 1, 0.5)
  )
  fit_dat$y <- with(fit_dat, 1 + x1 + 0.5 * a1 + rnorm(n_fit, sd = 0.2))

  new_dat <- data.frame(
    x1 = rnorm(n_new),
    a1 = rbinom(n_new, 1, 0.5)
  )

  out <- qstep(
    qmodel = qmodel,
    data = fit_dat,
    response = fit_dat["y"],
    newdata = new_dat,
    regime = rep(0L, n_new),
    txName = "a1"
  )

  expect_length(out$hats_mod, n_new)
  expect_length(out$hats_unmod, n_new)
})

test_that("qstep modifies predictions when regime differs from observed treatment", {
  qmodel <- modelObj::buildModelObj(
    model = ~ x1 + a1 + x1:a1,
    solver.method = "lm",
    predict.method = "predict.lm"
  )

  set.seed(3)
  n <- 150
  dat <- data.frame(
    x1 = rnorm(n),
    a1 = rep(c(0L, 1L), length.out = n)
  )
  dat$y <- with(dat, -0.2 + 0.9 * x1 + 1.5 * a1 + 0.6 * x1 * a1 + rnorm(n, sd = 0.1))

  regime <- rep(1L, n)
  out <- qstep(
    qmodel = qmodel,
    data = dat,
    response = dat["y"],
    newdata = dat,
    regime = regime,
    txName = "a1"
  )

  changed <- dat$a1 != regime
  expect_true(any(changed))
  expect_true(any(abs(out$hats_mod[changed] - out$hats_unmod[changed]) > 1e-8))
})
