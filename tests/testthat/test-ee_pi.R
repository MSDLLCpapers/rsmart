# ============================================================
# Tests for ee_pi.R: ee_psi_pi and ee_dpsi_pi
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# Helper -----------------------------------------------------------------------
make_pi_inputs <- function() {
  data("pcsttrial", package = "rsmart", envir = environment())
  df       <- pcsttrial
  df$t1    <- df$study_day_enroll
  df$t2    <- df$study_day_rerand
  df$t3    <- df$study_day_outcome
  K        <- 2L
  df$kappa <- get_kappa(df, max(df$t3), K)

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
  list(df = df, p_fits = piFitted[["p_fits"]], K = K)
}

# ee_psi_pi --------------------------------------------------------------------

test_that("ee_psi_pi returns a matrix with n rows", {
  inp <- make_pi_inputs()
  out <- ee_psi_pi(inp$df, inp$p_fits)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), nrow(inp$df))
})

test_that("ee_psi_pi has one column per propensity score parameter across stages", {
  inp      <- make_pi_inputs()
  out      <- ee_psi_pi(inp$df, inp$p_fits)
  # intercept-only models: 1 parameter per stage, 2 stages => 2 columns total
  n_params <- sum(vapply(inp$p_fits, function(f) length(f@fitObj$coefficients), integer(1L)))
  expect_equal(ncol(out), n_params)
})

test_that("ee_psi_pi column sums are zero at the MLE", {
  # Score equations equal zero at the MLE of a GLM
  inp  <- make_pi_inputs()
  out  <- ee_psi_pi(inp$df, inp$p_fits)
  sums <- colSums(out)
  expect_equal(sums, rep(0, ncol(out)), tolerance = 1e-8)
})

test_that("ee_psi_pi contains only finite values", {
  inp <- make_pi_inputs()
  out <- ee_psi_pi(inp$df, inp$p_fits)
  expect_true(all(is.finite(out)))
})

# ee_dpsi_pi -------------------------------------------------------------------

test_that("ee_dpsi_pi returns a square matrix with dimension equal to total pi parameters", {
  inp      <- make_pi_inputs()
  out      <- ee_dpsi_pi(inp$df, inp$p_fits)
  n_params <- sum(vapply(inp$p_fits, function(f) length(f@fitObj$coefficients), integer(1L)))
  expect_true(is.matrix(out) || inherits(out, "Matrix"))
  expect_equal(nrow(as.matrix(out)), n_params)
  expect_equal(ncol(as.matrix(out)), n_params)
})

test_that("ee_dpsi_pi is symmetric", {
  inp <- make_pi_inputs()
  out <- as.matrix(ee_dpsi_pi(inp$df, inp$p_fits))
  expect_equal(out, t(out), tolerance = 1e-10)
})

test_that("ee_dpsi_pi is negative semi-definite", {
  inp    <- make_pi_inputs()
  out    <- as.matrix(ee_dpsi_pi(inp$df, inp$p_fits))
  eigval <- eigen(out, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigval <= 1e-8))
})

test_that("ee_dpsi_pi contains only finite values", {
  inp <- make_pi_inputs()
  out <- ee_dpsi_pi(inp$df, inp$p_fits)
  expect_true(all(is.finite(as.vector(out))))
})

# ============================================================
# Additional tests with multi-covariate propensity models
# NOTE: gprime2 is defined but never called in both ee_psi_pi
# and ee_dpsi_pi, and ee_dpsi_pi contains a duplicate return
# statement. These are dead-code expressions that cannot be
# covered by any test — they would need to be removed from
# the source to reach 100% coverage.
# ============================================================

make_pi_inputs_cov <- function() {
  data("pcsttrial", package = "rsmart", envir = environment())
  df    <- pcsttrial
  df$t1 <- df$study_day_enroll
  df$t2 <- df$study_day_rerand
  df$t3 <- df$study_day_outcome
  K        <- 2L
  df$kappa <- get_kappa(df, max(df$t3), K)

  # Stage 1: intercept + one covariate (2 params)
  p1 <- modelObj::buildModelObj(
    model          = ~ comorbidity,
    solver.method  = "glm",
    solver.args    = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args   = list(type = "response")
  )
  # Stage 2: intercept + two covariates (3 params)
  p2 <- modelObj::buildModelObj(
    model          = ~ comorbidity + painmed,
    solver.method  = "glm",
    solver.args    = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args   = list(type = "response")
  )

  piFitted <- pi_fits(df, list(p1, p2))
  list(df = df, p_fits = piFitted[["p_fits"]], K = K)
}

test_that("ee_psi_pi works with multi-covariate models and has correct dimensions", {
  inp <- make_pi_inputs_cov()
  out <- ee_psi_pi(inp$df, inp$p_fits)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), nrow(inp$df))
  # 2 params (intercept + comorbidity) + 3 params (intercept + 2 covariates)
  n_params <- sum(vapply(inp$p_fits, function(f) length(f@fitObj$coefficients), integer(1L)))
  expect_equal(ncol(out), n_params)
})

test_that("ee_psi_pi with covariate models contains only finite values", {
  inp <- make_pi_inputs_cov()
  out <- ee_psi_pi(inp$df, inp$p_fits)
  expect_true(all(is.finite(out)))
})

test_that("ee_dpsi_pi with covariate models returns a square matrix of correct size", {
  inp      <- make_pi_inputs_cov()
  out      <- as.matrix(ee_dpsi_pi(inp$df, inp$p_fits))
  n_params <- sum(vapply(inp$p_fits, function(f) length(f@fitObj$coefficients), integer(1L)))
  expect_equal(dim(out), c(n_params, n_params))
})

test_that("ee_dpsi_pi with covariate models is negative semi-definite", {
  inp    <- make_pi_inputs_cov()
  out    <- as.matrix(ee_dpsi_pi(inp$df, inp$p_fits))
  eigval <- eigen(out, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigval <= 1e-8))
})

test_that("ee_dpsi_pi with covariate models contains only finite values", {
  inp <- make_pi_inputs_cov()
  out <- ee_dpsi_pi(inp$df, inp$p_fits)
  expect_true(all(is.finite(as.vector(out))))
})
