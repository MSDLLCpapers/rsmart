# ============================================================
# Tests for ee_value_derivatives_ipw.R:
#   ee_dpsiv_dpi_ipw, ee_dpsiv_dnu_ipw
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# Helper -----------------------------------------------------------------------
make_deriv_ipw_inputs <- function() {
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

  # Use regime 1 for per-regime derivative tests
  regime1 <- regime_all[[1]]

  list(
    df         = df,
    nus        = nus,
    pis        = pis,
    p_fits     = p_fits,
    regime_all = regime_all,
    regime1    = regime1,
    K          = K
  )
}

# ee_dpsiv_dpi_ipw -------------------------------------------------------------

test_that("ee_dpsiv_dpi_ipw returns a row vector (1 x p_k)", {
  inp <- make_deriv_ipw_inputs()
  out <- ee_dpsiv_dpi_ipw(
    df         = inp$df,
    pis        = inp$pis,
    nus        = inp$nus,
    regime_ind = inp$regime1$regime_ind,
    p_fits     = inp$p_fits,
    k          = 1L,
    qs         = NULL
  )
  n_params_k1 <- length(inp$p_fits[[1]]@fitObj$coefficients)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), 1L)
  expect_equal(ncol(out), n_params_k1)
})

test_that("ee_dpsiv_dpi_ipw contains finite values for stage 1", {
  inp <- make_deriv_ipw_inputs()
  out <- ee_dpsiv_dpi_ipw(
    df = inp$df, pis = inp$pis, nus = inp$nus,
    regime_ind = inp$regime1$regime_ind,
    p_fits = inp$p_fits, k = 1L, qs = NULL
  )
  expect_true(all(is.finite(out)))
})

test_that("ee_dpsiv_dpi_ipw contains finite values for stage 2", {
  inp <- make_deriv_ipw_inputs()
  out <- ee_dpsiv_dpi_ipw(
    df = inp$df, pis = inp$pis, nus = inp$nus,
    regime_ind = inp$regime1$regime_ind,
    p_fits = inp$p_fits, k = 2L, qs = NULL
  )
  expect_true(all(is.finite(out)))
})

test_that("ee_dpsiv_dpi_ipw differs across regimes", {
  inp  <- make_deriv_ipw_inputs()
  out1 <- ee_dpsiv_dpi_ipw(inp$df, inp$pis, inp$nus,
                           inp$regime_all[[1]]$regime_ind,
                           inp$p_fits, 1L, NULL)
  out4 <- ee_dpsiv_dpi_ipw(inp$df, inp$pis, inp$nus,
                           inp$regime_all[[4]]$regime_ind,
                           inp$p_fits, 1L, NULL)
  expect_false(isTRUE(all.equal(out1, out4)))
})

# ee_dpsiv_dnu_ipw -------------------------------------------------------------

test_that("ee_dpsiv_dnu_ipw returns a scalar for each k", {
  inp <- make_deriv_ipw_inputs()
  for (k in 1:(inp$K + 1L)) {
    out <- ee_dpsiv_dnu_ipw(
      df = inp$df, pis = inp$pis, nus = inp$nus,
      regime_ind = inp$regime1$regime_ind,
      k = k, qs = NULL
    )
    expect_length(out, 1L)
  }
})

test_that("ee_dpsiv_dnu_ipw is finite for all k", {
  inp <- make_deriv_ipw_inputs()
  for (k in 1:(inp$K + 1L)) {
    out <- ee_dpsiv_dnu_ipw(
      df = inp$df, pis = inp$pis, nus = inp$nus,
      regime_ind = inp$regime1$regime_ind,
      k = k, qs = NULL
    )
    expect_true(is.finite(out), label = paste("k =", k))
  }
})

test_that("ee_dpsiv_dnu_ipw is zero for k < K+1 (nu only affects IPW through nu_{K+1})", {
  inp <- make_deriv_ipw_inputs()
  for (k in 1:inp$K) {
    out <- ee_dpsiv_dnu_ipw(
      df = inp$df, pis = inp$pis, nus = inp$nus,
      regime_ind = inp$regime1$regime_ind,
      k = k, qs = NULL
    )
    expect_equal(out, 0, tolerance = 1e-12, label = paste("k =", k))
  }
})

test_that("ee_dpsiv_dnu_ipw is non-zero for k = K+1", {
  inp <- make_deriv_ipw_inputs()
  out <- ee_dpsiv_dnu_ipw(
    df = inp$df, pis = inp$pis, nus = inp$nus,
    regime_ind = inp$regime1$regime_ind,
    k = inp$K + 1L, qs = NULL
  )
  expect_false(out == 0)
})

test_that("ee_dpsiv_dnu_ipw result at k=K+1 is numeric", {
  inp <- make_deriv_ipw_inputs()
  out <- ee_dpsiv_dnu_ipw(
    df = inp$df, pis = inp$pis, nus = inp$nus,
    regime_ind = inp$regime1$regime_ind,
    k = inp$K + 1L, qs = NULL
  )
  expect_type(out, "double")
})
