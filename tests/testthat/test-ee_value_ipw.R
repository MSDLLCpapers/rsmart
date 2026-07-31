# ============================================================
# Tests for ee_value_ipw.R: ee_dpsiv_ipw
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# Helper -----------------------------------------------------------------------
make_value_ipw_inputs <- function() {
  data("pcsttrial", package = "rsmart", envir = environment())
  df       <- pcsttrial
  df$t1    <- df$study_day_enroll
  df$t2    <- df$study_day_rerand
  df$t3    <- df$study_day_outcome
  df$y     <- df$pctchange
  K        <- 2L
  df$kappa <- get_kappa(df, max(df$t3), K)
  nus      <- get_nu(df, K)

  p1 <- modelObj::buildModelObj(
    model          = ~ 1,
    solver.method  = "glm",
    solver.args    = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args   = list(type = "response")
  )
  p2 <- modelObj::buildModelObj(
    model          = ~ 1,
    solver.method  = "glm",
    solver.args    = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args   = list(type = "response")
  )
  piFitted <- pi_fits(df, list(p1, p2))
  pis      <- piFitted[["ps"]]
  p_fits   <- piFitted[["p_fits"]]

  regime_all <- regime_list_no_trt_resp(
    emb_regimes = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
    dat        = df,
    resp_trt    = list("r2" = list(0, 0, 0, 0))
  )

  # IPW: q_list = NULL, q_all = empty list
  vHats <- estimate_values(df, NULL, regime_all, TRUE, pis, nus)

  list(
    df         = df,
    nus        = nus,
    pis        = pis,
    p_fits     = p_fits,
    regime_all = regime_all,
    values     = vHats[["value"]],
    q_all      = list()
  )
}

# ee_dpsiv_ipw ----------------------------------------------------------------

test_that("ee_dpsiv_ipw returns a matrix", {
  inp <- make_value_ipw_inputs()
  out <- ee_dpsiv_ipw(inp$df, inp$pis, inp$nus, inp$regime_all, inp$p_fits,
                      inp$q_all, TRUE)
  expect_true(is.matrix(out))
})

test_that("ee_dpsiv_ipw has one row per regime", {
  inp <- make_value_ipw_inputs()
  out <- ee_dpsiv_ipw(inp$df, inp$pis, inp$nus, inp$regime_all, inp$p_fits,
                      inp$q_all, TRUE)
  expect_equal(nrow(out), length(inp$regime_all))
})

test_that("ee_dpsiv_ipw contains only finite values", {
  inp <- make_value_ipw_inputs()
  out <- ee_dpsiv_ipw(inp$df, inp$pis, inp$nus, inp$regime_all, inp$p_fits,
                      inp$q_all, TRUE)
  expect_true(all(is.finite(out)))
})

test_that("ee_dpsiv_ipw has fewer columns than ee_dpsiv (no beta params)", {
  # IPW has no Q-function parameters, so its column count should be smaller
  data("pcsttrial", package = "rsmart", envir = environment())
  df       <- pcsttrial
  df$t1    <- df$study_day_enroll
  df$t2    <- df$study_day_rerand
  df$t3    <- df$study_day_outcome
  df$y     <- df$pctchange
  K        <- 2L
  df$kappa <- get_kappa(df, max(df$t3), K)
  nus      <- get_nu(df, K)

  p1 <- modelObj::buildModelObj(
    model = ~ 1, solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )
  p2 <- modelObj::buildModelObj(
    model = ~ 1, solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )
  piFitted   <- pi_fits(df, list(p1, p2))
  pis        <- piFitted[["ps"]]
  p_fits     <- piFitted[["p_fits"]]
  regime_all <- regime_list_no_trt_resp(
    emb_regimes = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
    dat = df, resp_trt = list("r2" = list(0, 0, 0, 0))
  )

  q1 <- modelObj::buildModelObj(model = ~ a1, solver.method = "lm",
                                predict.method = "predict.lm")
  q2 <- modelObj::buildModelObj(model = ~ a2, solver.method = "lm",
                                predict.method = "predict.lm")
  q_list <- list(q1, q2)

  vHats_aug <- estimate_values(df, q_list, regime_all, TRUE, pis, nus)
  vHats_ipw <- estimate_values(df, NULL,   regime_all, TRUE, pis, nus)

  out_aug <- ee_dpsiv(df, pis, nus, regime_all, p_fits, vHats_aug[["q_all"]],
                      TRUE, q_list)
  out_ipw <- ee_dpsiv_ipw(df, pis, nus, regime_all, p_fits, list(), TRUE)

  expect_lt(ncol(out_ipw), ncol(out_aug))
})

test_that("ee_dpsiv_ipw last L columns form a negative diagonal block", {
  inp <- make_value_ipw_inputs()
  out <- ee_dpsiv_ipw(inp$df, inp$pis, inp$nus, inp$regime_all, inp$p_fits,
                      inp$q_all, TRUE)
  L        <- length(inp$values)
  last_blk <- out[, (ncol(out) - L + 1):ncol(out)]
  expected <- -inp$nus$ns * diag(L)
  expect_equal(last_blk, expected, tolerance = 1e-10)
})
