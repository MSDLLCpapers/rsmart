# ============================================================
# Tests for ee_beta.R: ee_psi_beta and ee_dpsi_beta
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# Helpers ----------------------------------------------------------------------
# These functions build the expensive pipeline (data → kappa/nu → pi_fits →
# regime list → estimate_values). They are called exactly once each at file
# load time (see "Shared fixtures" section below), so no test re-runs the
# pipeline on its own.

make_beta_inputs <- function() {
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

  list(
    df         = df,
    nus        = nus,
    regime_all = regime_all,
    q_list     = q_list,
    q_all      = vHats[["q_all"]]
  )
}

# Helper: feasible_sets_indicator = FALSE (all stages fully randomised) ----------
make_beta_inputs_nofsi <- function() {
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
  pis <- pi_fits(df, list(p1, p2))[["ps"]]

  regime_all <- regime_list_no_trt_resp(
    emb_regimes = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
    dat = df, resp_trt = list("r2" = list(0, 0, 0, 0))
  )

  q1 <- modelObj::buildModelObj(model = ~ a1, solver.method = "lm",
                                predict.method = "predict.lm")
  q2 <- modelObj::buildModelObj(model = ~ a2, solver.method = "lm",
                                predict.method = "predict.lm")
  q_list <- list(q1, q2)

  # Pass FALSE so get_q_fits uses the no-feasible-sets code path
  vHats <- estimate_values(df, q_list, regime_all, FALSE, pis, nus)

  list(df = df, nus = nus, regime_all = regime_all,
       q_list = q_list, q_all = vHats[["q_all"]])
}

# Helper: split stage-2 models (length(q_list[[2]]) > 1) ----------------------
make_beta_inputs_split <- function() {
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
  pis <- pi_fits(df, list(p1, p2))[["ps"]]

  regime_all <- regime_list_no_trt_resp(
    emb_regimes = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
    dat = df, resp_trt = list("r2" = list(0, 0, 0, 0))
  )

  q1     <- modelObj::buildModelObj(model = ~ a1, solver.method = "lm",
                                    predict.method = "predict.lm")
  q2_r0  <- modelObj::buildModelObj(model = ~ a2, solver.method = "lm",
                                    predict.method = "predict.lm")
  q2_r1  <- modelObj::buildModelObj(model = ~ 1,  solver.method = "lm",
                                    predict.method = "predict.lm")
  q_list <- list(q1, list(r0 = q2_r0, r1 = q2_r1))

  vHats <- estimate_values(df, q_list, regime_all, TRUE, pis, nus)

  list(df = df, nus = nus, regime_all = regime_all,
       q_list = q_list, q_all = vHats[["q_all"]])
}

# Shared fixtures: each pipeline is run exactly once per file load -------------
# All test_that() blocks reference these objects; no helper is called inside a
# test.  If any fixture fails to build, only that one construction fails with a
# clear error before any test runs.
.inp       <- make_beta_inputs()
.inp_nofsi <- make_beta_inputs_nofsi()
.inp_split <- make_beta_inputs_split()

# ee_psi_beta ------------------------------------------------------------------

test_that("ee_psi_beta returns a matrix with n rows", {
  out <- ee_psi_beta(.inp$df, .inp$q_all, TRUE, .inp$q_list)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), nrow(.inp$df))
})

test_that("ee_psi_beta has one column per beta parameter across all regimes and stages", {
  out <- ee_psi_beta(.inp$df, .inp$q_all, TRUE, .inp$q_list)
  # Count actual fitted coefficients from the first regime's q_fits
  # (all regimes share the same model structure)
  n_params_per_regime <- sum(vapply(
    .inp$q_all[[1]]$q_fits,
    function(f) length(f@fitObj$coefficients),
    integer(1L)
  ))
  expect_equal(ncol(out), length(.inp$q_all) * n_params_per_regime)
})

test_that("ee_psi_beta column sums are zero at the OLS estimates", {
  # Score equations of OLS are zero at the MLE by the normal equations
  out  <- ee_psi_beta(.inp$df, .inp$q_all, TRUE, .inp$q_list)
  sums <- colSums(out)
  expect_equal(unname(sums), rep(0, ncol(out)), tolerance = 1e-8)
})

test_that("ee_psi_beta contains only finite values", {
  out <- ee_psi_beta(.inp$df, .inp$q_all, TRUE, .inp$q_list)
  expect_true(all(is.finite(out)))
})

# ee_dpsi_beta -----------------------------------------------------------------

test_that("ee_dpsi_beta returns a square matrix", {
  out <- as.matrix(ee_dpsi_beta(.inp$df, .inp$q_all, .inp$regime_all, TRUE, .inp$q_list))
  expect_equal(nrow(out), ncol(out))
})

test_that("ee_dpsi_beta dimension matches total beta parameter count", {
  out      <- as.matrix(ee_dpsi_beta(.inp$df, .inp$q_all, .inp$regime_all, TRUE, .inp$q_list))
  psi_beta <- ee_psi_beta(.inp$df, .inp$q_all, TRUE, .inp$q_list)
  expect_equal(nrow(out), ncol(psi_beta))
})

test_that("ee_dpsi_beta is negative semi-definite", {
  # Each diagonal block is -X^T X, which is negative semi-definite
  out    <- as.matrix(ee_dpsi_beta(.inp$df, .inp$q_all, .inp$regime_all, TRUE, .inp$q_list))
  eigval <- eigen(out, symmetric = FALSE, only.values = TRUE)$values
  expect_true(all(Re(eigval) <= 1e-6))
})

test_that("ee_dpsi_beta contains only finite values", {
  out <- ee_dpsi_beta(.inp$df, .inp$q_all, .inp$regime_all, TRUE, .inp$q_list)
  expect_true(all(is.finite(as.vector(out))))
})

# ee_psi_beta: feasible_sets_indicator = FALSE ------------------------------------

test_that("ee_psi_beta returns a matrix when feasible_sets_indicator is FALSE", {
  out <- ee_psi_beta(.inp_nofsi$df, .inp_nofsi$q_all, FALSE, .inp_nofsi$q_list)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), nrow(.inp_nofsi$df))
})

test_that("ee_psi_beta column sums are zero when feasible_sets_indicator is FALSE", {
  out  <- ee_psi_beta(.inp_nofsi$df, .inp_nofsi$q_all, FALSE, .inp_nofsi$q_list)
  expect_equal(unname(colSums(out)), rep(0, ncol(out)), tolerance = 1e-8)
})

test_that("ee_psi_beta contains only finite values when feasible_sets_indicator is FALSE", {
  out <- ee_psi_beta(.inp_nofsi$df, .inp_nofsi$q_all, FALSE, .inp_nofsi$q_list)
  expect_true(all(is.finite(out)))
})

# ee_psi_beta: split stage-2 models ---------------------------------------------

test_that("ee_psi_beta returns a matrix with split stage-2 models", {
  out <- ee_psi_beta(.inp_split$df, .inp_split$q_all, TRUE, .inp_split$q_list)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), nrow(.inp_split$df))
})

test_that("ee_psi_beta with split models has more columns than single-model version", {
  out_single <- ee_psi_beta(.inp$df,       .inp$q_all,       TRUE, .inp$q_list)
  out_split  <- ee_psi_beta(.inp_split$df, .inp_split$q_all, TRUE, .inp_split$q_list)
  expect_gt(ncol(out_split), ncol(out_single))
})

test_that("ee_psi_beta with split models contains only finite values", {
  out <- ee_psi_beta(.inp_split$df, .inp_split$q_all, TRUE, .inp_split$q_list)
  expect_true(all(is.finite(out)))
})

# ee_dpsi_beta: feasible_sets_indicator = FALSE -----------------------------------

test_that("ee_dpsi_beta returns a square matrix when feasible_sets_indicator is FALSE", {
  out <- as.matrix(ee_dpsi_beta(.inp_nofsi$df, .inp_nofsi$q_all, .inp_nofsi$regime_all, FALSE, .inp_nofsi$q_list))
  expect_equal(nrow(out), ncol(out))
})

test_that("ee_dpsi_beta is negative semi-definite when feasible_sets_indicator is FALSE", {
  out    <- as.matrix(ee_dpsi_beta(.inp_nofsi$df, .inp_nofsi$q_all, .inp_nofsi$regime_all, FALSE, .inp_nofsi$q_list))
  eigval <- eigen(out, symmetric = FALSE, only.values = TRUE)$values
  expect_true(all(Re(eigval) <= 1e-6))
})

test_that("ee_dpsi_beta contains only finite values when feasible_sets_indicator is FALSE", {
  out <- ee_dpsi_beta(.inp_nofsi$df, .inp_nofsi$q_all, .inp_nofsi$regime_all, FALSE, .inp_nofsi$q_list)
  expect_true(all(is.finite(as.vector(out))))
})

# ee_dpsi_beta: split stage-2 models --------------------------------------------

test_that("ee_dpsi_beta returns a square matrix with split stage-2 models", {
  out <- as.matrix(ee_dpsi_beta(.inp_split$df, .inp_split$q_all, .inp_split$regime_all, TRUE, .inp_split$q_list))
  expect_equal(nrow(out), ncol(out))
})

test_that("ee_dpsi_beta with split models contains only finite values", {
  out <- ee_dpsi_beta(.inp_split$df, .inp_split$q_all, .inp_split$regime_all, TRUE, .inp_split$q_list)
  expect_true(all(is.finite(as.vector(out))))
})

test_that("ee_dpsi_beta with split models is larger than the single-model version", {
  out_single <- as.matrix(ee_dpsi_beta(.inp$df,       .inp$q_all,       .inp$regime_all,       TRUE, .inp$q_list))
  out_split  <- as.matrix(ee_dpsi_beta(.inp_split$df, .inp_split$q_all, .inp_split$regime_all, TRUE, .inp_split$q_list))
  expect_gt(nrow(out_split), nrow(out_single))
})
