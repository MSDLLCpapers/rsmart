# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

make_getQfits_data <- function(n = 200, with_r2 = FALSE) {
  set.seed(123)
  dat <- data.frame(
    x11 = rnorm(n),
    x12 = rnorm(n),
    x21 = rnorm(n),
    a1 = rbinom(n, 1, 0.5),
    a2 = rbinom(n, 1, 0.5)
  )

  dat$y <- with(dat, 0.5 + x11 - 0.4 * x12 + 0.3 * x21 + 0.2 * a1 - 0.1 * a2 + rnorm(n, sd = 0.2))

  if (with_r2) {
    dat$r2 <- rbinom(n, 1, 0.5)
    dat$kappa <- 3L
    # ensure responder stage-2 prediction branch (new = kappa == 2) has data
    idx_resp <- which(dat$r2 == 1)
    dat$kappa[idx_resp[seq_len(min(20, length(idx_resp)))]] <- 2L
  } else {
    dat$kappa <- 3L
  }

  dat
}

make_q_models <- function(split_stage2 = FALSE) {
  q1 <- modelObj::buildModelObj(
    model = ~ x11 + x12 + a1 + x11:a1 + x12:a1,
    solver.method = "lm",
    predict.method = "predict.lm"
  )

  q2 <- modelObj::buildModelObj(
    model = ~ x11 + x12 + x21 + a1 + a2 + a1:a2 + x11:a1 + x12:a1,
    solver.method = "lm",
    predict.method = "predict.lm"
  )

  if (!split_stage2) {
    return(list(q1, q2))
  }

  q2r <- modelObj::buildModelObj(
    model = ~ x11 + x12 + x21 + a1 + x11:a1 + x12:a1,
    solver.method = "lm",
    predict.method = "predict.lm"
  )

  list(q1, list(r0 = q2, r1 = q2r))
}

test_that("get_q_fits works for no feasible sets", {
  dat <- make_getQfits_data(with_r2 = FALSE)
  q_list <- make_q_models(split_stage2 = FALSE)
  regime <- as.matrix(dat[, c("a1", "a2")])

  out <- get_q_fits(dat, q_list, regime, feasible_sets_indicator = FALSE)

  expect_named(out, c("q_fits", "mod_regime_vhats", "unmod_regime_vhats"))
  expect_length(out$q_fits, 2)
  expect_equal(dim(out$mod_regime_vhats), c(nrow(dat), 3))
  expect_equal(dim(out$unmod_regime_vhats), c(nrow(dat), 3))
  expect_equal(colnames(out$mod_regime_vhats), c("q1", "q2", "q3"))
  expect_equal(colnames(out$unmod_regime_vhats), c("q1_nochange", "q2_nochange", "q3_nochange"))
})

test_that("get_q_fits works for feasible sets without response-model split", {
  dat <- make_getQfits_data(with_r2 = FALSE)
  q_list <- make_q_models(split_stage2 = FALSE)
  regime <- as.matrix(dat[, c("a1", "a2")])

  out <- get_q_fits(dat, q_list, regime, feasible_sets_indicator = TRUE)

  expect_named(out, c("q_fits", "mod_regime_vhats", "unmod_regime_vhats"))
  expect_length(out$q_fits, 2)
  expect_false(is.list(out$q_fits[[2]]) && all(c("r0", "r1") %in% names(out$q_fits[[2]])))
  expect_true(all(is.finite(out$mod_regime_vhats)))
  expect_true(all(is.finite(out$unmod_regime_vhats)))
})

test_that("get_q_fits works for feasible sets with stage-2 r0/r1 split models", {
  dat <- make_getQfits_data(with_r2 = TRUE)
  q_list <- make_q_models(split_stage2 = TRUE)
  regime <- as.matrix(dat[, c("a1", "a2")])

  out <- get_q_fits(dat, q_list, regime, feasible_sets_indicator = TRUE)

  expect_named(out, c("q_fits", "mod_regime_vhats", "unmod_regime_vhats"))
  expect_length(out$q_fits, 2)
  expect_true(is.list(out$q_fits[[2]]))
  expect_named(out$q_fits[[2]], c("r0", "r1"))
  expect_true(all(is.finite(out$mod_regime_vhats)))
  expect_true(all(is.finite(out$unmod_regime_vhats)))
})
