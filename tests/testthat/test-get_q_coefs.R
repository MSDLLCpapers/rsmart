# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

make_lm_fit <- function(y, x = NULL) {
  dat <- data.frame(y = y)
  if (is.null(x)) {
    mod <- modelObj::buildModelObj(
      model = ~ 1,
      solver.method = "lm",
      predict.method = "predict.lm"
    )
  } else {
    dat$x <- x
    mod <- modelObj::buildModelObj(
      model = ~ x,
      solver.method = "lm",
      predict.method = "predict.lm"
    )
  }
  modelObj::fit(object = mod, data = dat, response = dat["y"])
}

test_that("get_q_coefs concatenates coefficients for no feasible sets structure", {
  fit_k1 <- make_lm_fit(y = 1:10, x = 0:9)
  fit_k2 <- make_lm_fit(y = 2:11, x = 1:10)

  q_all <- list(
    list(q_fits = list(fit_k1, fit_k2))
  )

  expected <- c(fit_k1@fitObj$coefficients,
                fit_k2@fitObj$coefficients)

  expect_equal(get_q_coefs(q_all), expected)
})

test_that("get_q_coefs concatenates coefficients for feasible sets AIPW structure", {
  fit_reg1_k1 <- make_lm_fit(y = 3:12, x = 1:10)
  fit_reg1_k2 <- make_lm_fit(y = 5:14, x = 0:9)
  fit_reg2_k1 <- make_lm_fit(y = 2:11, x = 2:11)

  q_all <- list(
    list(q_fits = list(fit_reg1_k1, fit_reg1_k2)),
    list(q_fits = list(fit_reg2_k1))
  )

  expected <- c(fit_reg1_k1@fitObj$coefficients,
                fit_reg1_k2@fitObj$coefficients,
                fit_reg2_k1@fitObj$coefficients)

  expect_equal(get_q_coefs(q_all), expected)
})

test_that("get_q_coefs concatenates coefficients for feasible sets IAIPW split models", {
  fit_k1 <- make_lm_fit(y = 1:10, x = 0:9)
  fit_k2_r0 <- make_lm_fit(y = 4:13, x = 1:10)
  fit_k2_r1 <- make_lm_fit(y = 7:16, x = 1:10)

  q_all <- list(
    list(q_fits = list(fit_k1, list(r0 = fit_k2_r0, r1 = fit_k2_r1)))
  )

  expected <- c(fit_k1@fitObj$coefficients,
                unlist(lapply(list(r0 = fit_k2_r0, r1 = fit_k2_r1),
                              function(x) x@fitObj$coefficients)))

  expect_equal(get_q_coefs(q_all), expected)
})
