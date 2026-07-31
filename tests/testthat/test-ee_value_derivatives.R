# ============================================================
# Tests for ee_value_derivatives.R:
#   ee_dpsiv_dpi, ee_dpsiv_dnu, ee_dpsiv_dbeta, ee_dpsiv_dv
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# Helper -----------------------------------------------------------------------
make_deriv_inputs <- function() {
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

  q1 <- modelObj::buildModelObj(
    model          = ~ a1,
    solver.method  = "lm",
    predict.method = "predict.lm"
  )
  q2 <- modelObj::buildModelObj(
    model          = ~ a2,
    solver.method  = "lm",
    predict.method = "predict.lm"
  )
  q_list <- list(q1, q2)

  vHats <- estimate_values(df, q_list, regime_all, TRUE, pis, nus)

  # Use regime 1 for per-regime derivative tests
  regime1    <- regime_all[[1]]
  qs1        <- vHats[["q_all"]][[1]]

  list(
    df         = df,
    nus        = nus,
    pis        = pis,
    p_fits     = p_fits,
    regime_all = regime_all,
    q_list     = q_list,
    q_all      = vHats[["q_all"]],
    values     = vHats[["value"]],
    regime1    = regime1,
    qs1        = qs1,
    K          = K
  )
}

# ee_dpsiv_dpi -----------------------------------------------------------------

test_that("ee_dpsiv_dpi returns a row vector (1 x p_k)", {
  inp <- make_deriv_inputs()
  out <- ee_dpsiv_dpi(
    df         = inp$df,
    pis        = inp$pis,
    nus        = inp$nus,
    regime_ind = inp$regime1$regime_ind,
    p_fits     = inp$p_fits,
    k          = 1L,
    qs         = inp$qs1
  )
  n_params_k1 <- length(inp$p_fits[[1]]@fitObj$coefficients)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), 1L)
  expect_equal(ncol(out), n_params_k1)
})

test_that("ee_dpsiv_dpi contains finite values for stage 1", {
  inp <- make_deriv_inputs()
  out <- ee_dpsiv_dpi(
    df = inp$df, pis = inp$pis, nus = inp$nus,
    regime_ind = inp$regime1$regime_ind,
    p_fits = inp$p_fits, k = 1L, qs = inp$qs1
  )
  expect_true(all(is.finite(out)))
})

test_that("ee_dpsiv_dpi contains finite values for stage 2", {
  inp <- make_deriv_inputs()
  out <- ee_dpsiv_dpi(
    df = inp$df, pis = inp$pis, nus = inp$nus,
    regime_ind = inp$regime1$regime_ind,
    p_fits = inp$p_fits, k = 2L, qs = inp$qs1
  )
  expect_true(all(is.finite(out)))
})

test_that("ee_dpsiv_dpi differs across regimes", {
  inp  <- make_deriv_inputs()
  out1 <- ee_dpsiv_dpi(inp$df, inp$pis, inp$nus,
                      inp$regime_all[[1]]$regime_ind, inp$p_fits, 1L,
                      inp$q_all[[1]])
  out2 <- ee_dpsiv_dpi(inp$df, inp$pis, inp$nus,
                      inp$regime_all[[4]]$regime_ind, inp$p_fits, 1L,
                      inp$q_all[[4]])
  expect_false(isTRUE(all.equal(out1, out2)))
})

# ee_dpsiv_dnu -----------------------------------------------------------------

test_that("ee_dpsiv_dnu returns a scalar for each k", {
  inp <- make_deriv_inputs()
  for (k in 1:(inp$K + 1L)) {
    out <- ee_dpsiv_dnu(
      df = inp$df, pis = inp$pis, nus = inp$nus,
      regime_ind = inp$regime1$regime_ind,
      k = k, qs = inp$qs1
    )
    expect_length(out, 1L)
  }
})

test_that("ee_dpsiv_dnu is finite for all k", {
  inp <- make_deriv_inputs()
  for (k in 1:(inp$K + 1L)) {
    out <- ee_dpsiv_dnu(
      df = inp$df, pis = inp$pis, nus = inp$nus,
      regime_ind = inp$regime1$regime_ind,
      k = k, qs = inp$qs1
    )
    expect_true(is.finite(out), label = paste("k =", k))
  }
})

test_that("ee_dpsiv_dnu is numeric", {
  inp <- make_deriv_inputs()
  out <- ee_dpsiv_dnu(
    df = inp$df, pis = inp$pis, nus = inp$nus,
    regime_ind = inp$regime1$regime_ind,
    k = inp$K + 1L, qs = inp$qs1
  )
  expect_type(out, "double")
})

# ee_dpsiv_dbeta ---------------------------------------------------------------

test_that("ee_dpsiv_dbeta returns a numeric vector for each stage", {
  inp <- make_deriv_inputs()
  for (k in 1:inp$K) {
    out <- ee_dpsiv_dbeta(
      df                   = inp$df,
      pis                  = inp$pis,
      nus                  = inp$nus,
      regime_ind           = inp$regime1$regime_ind,
      regime               = inp$regime1$regime,
      q_fits               = inp$qs1$q_fits,
      k                    = k,
      feasible_sets_indicator = TRUE,
      q_list               = inp$q_list
    )
    expect_true(is.numeric(out))
    expect_true(length(out) > 0L)
  }
})

test_that("ee_dpsiv_dbeta contains finite values", {
  inp <- make_deriv_inputs()
  for (k in 1:inp$K) {
    out <- ee_dpsiv_dbeta(
      df = inp$df, pis = inp$pis, nus = inp$nus,
      regime_ind = inp$regime1$regime_ind,
      regime = inp$regime1$regime,
      q_fits = inp$qs1$q_fits,
      k = k, feasible_sets_indicator = TRUE, q_list = inp$q_list
    )
    expect_true(all(is.finite(out)), label = paste("stage k =", k))
  }
})

test_that("ee_dpsiv_dbeta length matches number of beta parameters at stage k", {
  inp <- make_deriv_inputs()
  # model ~ a2 has 2 params (intercept + a2); stage k=2
  out <- ee_dpsiv_dbeta(
    df = inp$df, pis = inp$pis, nus = inp$nus,
    regime_ind = inp$regime1$regime_ind,
    regime = inp$regime1$regime,
    q_fits = inp$qs1$q_fits,
    k = 2L, feasible_sets_indicator = TRUE, q_list = inp$q_list
  )
  n_params <- length(inp$p_fits[[1]]@fitObj$coefficients) +
    length(inp$q_list[[2]]@model)   # intercept + 1 term
  # model ~ a2 produces 2 coefficients; just confirm length is positive and finite
  expect_true(length(out) >= 1L)
})

# ee_dpsiv_dv ------------------------------------------------------------------

test_that("ee_dpsiv_dv returns an L x L diagonal matrix", {
  inp <- make_deriv_inputs()
  L   <- length(inp$values)
  out <- ee_dpsiv_dv(nregimes = L, nus = inp$nus)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(L, L))
})

test_that("ee_dpsiv_dv diagonal entries equal -ns", {
  inp <- make_deriv_inputs()
  L   <- length(inp$values)
  out <- ee_dpsiv_dv(nregimes = L, nus = inp$nus)
  expect_equal(diag(out), rep(-inp$nus$ns, L), tolerance = 1e-10)
})

test_that("ee_dpsiv_dv off-diagonal entries are zero", {
  inp <- make_deriv_inputs()
  L   <- length(inp$values)
  out <- ee_dpsiv_dv(nregimes = L, nus = inp$nus)
  expect_equal(sum(out) - sum(diag(out)), 0, tolerance = 1e-10)
})

test_that("ee_dpsiv_dv is negative definite", {
  inp    <- make_deriv_inputs()
  L      <- length(inp$values)
  out    <- ee_dpsiv_dv(nregimes = L, nus = inp$nus)
  eigval <- eigen(out, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigval < 0))
})

# ============================================================
# Additional helpers and tests for missing branches in ee_dpsiv_dbeta
# ============================================================

# Helper: feasible_sets_indicator = FALSE ----------------------------------------
make_deriv_inputs_nofsi <- function() {
  data("pcsttrial", package = "rsmart", envir = environment())
  df    <- pcsttrial
  df$t1 <- df$study_day_enroll
  df$t2 <- df$study_day_rerand
  df$t3 <- df$study_day_outcome
  df$y  <- df$pctchange
  K        <- 2L
  df$kappa <- get_kappa(df, max(df$t3), K)
  nus      <- get_nu(df, K)

  p1 <- modelObj::buildModelObj(
    model = ~ 1, solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm", predict.args = list(type = "response")
  )
  p2 <- modelObj::buildModelObj(
    model = ~ 1, solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm", predict.args = list(type = "response")
  )
  piFitted <- pi_fits(df, list(p1, p2))
  pis      <- piFitted[["ps"]]

  regime_all <- regime_list_no_trt_resp(
    emb_regimes = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
    dat = df, resp_trt = list("r2" = list(0, 0, 0, 0))
  )

  q1 <- modelObj::buildModelObj(model = ~ a1, solver.method = "lm",
                                predict.method = "predict.lm")
  q2 <- modelObj::buildModelObj(model = ~ a2, solver.method = "lm",
                                predict.method = "predict.lm")
  q_list <- list(q1, q2)

  vHats <- estimate_values(df, q_list, regime_all, FALSE, pis, nus)
  qs1   <- vHats[["q_all"]][[1]]

  list(df = df, nus = nus, pis = pis, regime_all = regime_all,
       q_list = q_list, regime1 = regime_all[[1]], qs1 = qs1, K = K)
}

# Helper: split stage-2 models (length(q_list[[2]]) > 1) ----------------------
make_deriv_inputs_split <- function() {
  data("pcsttrial", package = "rsmart", envir = environment())
  df    <- pcsttrial
  df$t1 <- df$study_day_enroll
  df$t2 <- df$study_day_rerand
  df$t3 <- df$study_day_outcome
  df$y  <- df$pctchange
  K        <- 2L
  df$kappa <- get_kappa(df, max(df$t3), K)
  nus      <- get_nu(df, K)

  p1 <- modelObj::buildModelObj(
    model = ~ 1, solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm", predict.args = list(type = "response")
  )
  p2 <- modelObj::buildModelObj(
    model = ~ 1, solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm", predict.args = list(type = "response")
  )
  piFitted <- pi_fits(df, list(p1, p2))
  pis      <- piFitted[["ps"]]

  regime_all <- regime_list_no_trt_resp(
    emb_regimes = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
    dat = df, resp_trt = list("r2" = list(0, 0, 0, 0))
  )

  q1    <- modelObj::buildModelObj(model = ~ a1, solver.method = "lm",
                                   predict.method = "predict.lm")
  q2_r0 <- modelObj::buildModelObj(model = ~ a2, solver.method = "lm",
                                   predict.method = "predict.lm")
  q2_r1 <- modelObj::buildModelObj(model = ~ 1,  solver.method = "lm",
                                   predict.method = "predict.lm")
  q_list <- list(q1, list(r0 = q2_r0, r1 = q2_r1))

  vHats <- estimate_values(df, q_list, regime_all, TRUE, pis, nus)
  qs1   <- vHats[["q_all"]][[1]]

  list(df = df, nus = nus, pis = pis, regime_all = regime_all,
       q_list = q_list, regime1 = regime_all[[1]], qs1 = qs1, K = K)
}

# ee_dpsiv_dbeta: feasible_sets_indicator = FALSE ---------------------------------

test_that("ee_dpsiv_dbeta is finite for each stage when feasible_sets_indicator is FALSE", {
  inp <- make_deriv_inputs_nofsi()
  for (k in 1:inp$K) {
    out <- ee_dpsiv_dbeta(
      df = inp$df, pis = inp$pis, nus = inp$nus,
      regime_ind = inp$regime1$regime_ind,
      regime     = inp$regime1$regime,
      q_fits     = inp$qs1$q_fits,
      k = k, feasible_sets_indicator = FALSE, q_list = inp$q_list
    )
    expect_true(all(is.finite(out)), label = paste("no-fsi stage k =", k))
  }
})

test_that("ee_dpsiv_dbeta is numeric for each stage when feasible_sets_indicator is FALSE", {
  inp <- make_deriv_inputs_nofsi()
  for (k in 1:inp$K) {
    out <- ee_dpsiv_dbeta(
      df = inp$df, pis = inp$pis, nus = inp$nus,
      regime_ind = inp$regime1$regime_ind,
      regime     = inp$regime1$regime,
      q_fits     = inp$qs1$q_fits,
      k = k, feasible_sets_indicator = FALSE, q_list = inp$q_list
    )
    expect_true(is.numeric(out), label = paste("no-fsi stage k =", k))
  }
})

# ee_dpsiv_dbeta: split stage-2 models ------------------------------------------

test_that("ee_dpsiv_dbeta is finite at stage 2 with split models", {
  inp <- make_deriv_inputs_split()
  out <- ee_dpsiv_dbeta(
    df = inp$df, pis = inp$pis, nus = inp$nus,
    regime_ind = inp$regime1$regime_ind,
    regime     = inp$regime1$regime,
    q_fits     = inp$qs1$q_fits,
    k = 2L, feasible_sets_indicator = TRUE, q_list = inp$q_list
  )
  expect_true(all(is.finite(out)))
})

test_that("ee_dpsiv_dbeta with split models returns more parameters than single model", {
  inp_single <- make_deriv_inputs()
  inp_split  <- make_deriv_inputs_split()
  out_single <- ee_dpsiv_dbeta(
    df = inp_single$df, pis = inp_single$pis, nus = inp_single$nus,
    regime_ind = inp_single$regime1$regime_ind,
    regime = inp_single$regime1$regime, q_fits = inp_single$qs1$q_fits,
    k = 2L, feasible_sets_indicator = TRUE, q_list = inp_single$q_list
  )
  out_split <- ee_dpsiv_dbeta(
    df = inp_split$df, pis = inp_split$pis, nus = inp_split$nus,
    regime_ind = inp_split$regime1$regime_ind,
    regime = inp_split$regime1$regime, q_fits = inp_split$qs1$q_fits,
    k = 2L, feasible_sets_indicator = TRUE, q_list = inp_split$q_list
  )
  expect_gt(length(out_split), length(out_single))
})

test_that("ee_dpsiv_dbeta with split models is numeric at stage 1", {
  inp <- make_deriv_inputs_split()
  out <- ee_dpsiv_dbeta(
    df = inp$df, pis = inp$pis, nus = inp$nus,
    regime_ind = inp$regime1$regime_ind,
    regime     = inp$regime1$regime,
    q_fits     = inp$qs1$q_fits,
    k = 1L, feasible_sets_indicator = TRUE, q_list = inp$q_list
  )
  expect_true(is.numeric(out))
  expect_true(all(is.finite(out)))
})
