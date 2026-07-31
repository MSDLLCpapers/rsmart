# ============================================================
# Tests for ee_value.R: ee_psi_v and ee_dpsiv
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# Helper -----------------------------------------------------------------------
make_value_inputs <- function() {
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

  vHats  <- estimate_values(df, q_list, regime_all, TRUE, pis, nus)

  list(
    df         = df,
    nus        = nus,
    pis        = pis,
    p_fits     = p_fits,
    regime_all = regime_all,
    q_list     = q_list,
    values     = vHats[["value"]],
    dfs        = vHats[["df"]],
    q_all      = vHats[["q_all"]]
  )
}

# ee_psi_v ---------------------------------------------------------------------

test_that("ee_psi_v returns a matrix with n rows and L columns", {
  inp   <- make_value_inputs()
  cKap  <- (inp$df$kappa > 0)
  out   <- ee_psi_v(inp$dfs, inp$values, cKap)
  L     <- length(inp$values)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), nrow(inp$df))
  expect_equal(ncol(out), L)
})

test_that("ee_psi_v column sums are zero at estimated values", {
  inp  <- make_value_inputs()
  cKap <- (inp$df$kappa > 0)
  out  <- ee_psi_v(inp$dfs, inp$values, cKap)
  sums <- colSums(out)
  expect_equal(sums, rep(0, ncol(out)), tolerance = 1e-8)
})

test_that("ee_psi_v contains only finite values", {
  inp  <- make_value_inputs()
  cKap <- (inp$df$kappa > 0)
  out  <- ee_psi_v(inp$dfs, inp$values, cKap)
  expect_true(all(is.finite(out)))
})

test_that("ee_psi_v column for regime ell sums to zero when value is mean of terms", {
  # Value is defined as sum(vTerms) / ns, so sum(vTerms) - value * ns = 0
  inp  <- make_value_inputs()
  cKap <- (inp$df$kappa > 0)
  out  <- ee_psi_v(inp$dfs, inp$values, cKap)
  for (ell in seq_along(inp$values)) {
    expect_equal(sum(out[, ell]), 0, tolerance = 1e-8,
                 label = paste("Column sum for regime", ell))
  }
})

# ee_dpsiv --------------------------------------------------------------------

test_that("ee_dpsiv returns a matrix with L rows", {
  inp <- make_value_inputs()
  out <- ee_dpsiv(inp$df, inp$pis, inp$nus, inp$regime_all, inp$p_fits,
                  inp$q_all, TRUE, inp$q_list)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), length(inp$values))
})

test_that("ee_dpsiv row count matches number of regimes", {
  inp <- make_value_inputs()
  out <- ee_dpsiv(inp$df, inp$pis, inp$nus, inp$regime_all, inp$p_fits,
                  inp$q_all, TRUE, inp$q_list)
  expect_equal(nrow(out), length(inp$regime_all))
})

test_that("ee_dpsiv contains only finite values", {
  inp <- make_value_inputs()
  out <- ee_dpsiv(inp$df, inp$pis, inp$nus, inp$regime_all, inp$p_fits,
                  inp$q_all, TRUE, inp$q_list)
  expect_true(all(is.finite(out)))
})

test_that("ee_dpsiv last L columns form a negative diagonal block", {
  # The last L columns correspond to -dPsiV/dV = -ns * I_L
  inp <- make_value_inputs()
  out <- ee_dpsiv(inp$df, inp$pis, inp$nus, inp$regime_all, inp$p_fits,
                  inp$q_all, TRUE, inp$q_list)
  L        <- length(inp$values)
  last_blk <- out[, (ncol(out) - L + 1):ncol(out)]
  expected <- -inp$nus$ns * diag(L)
  expect_equal(last_blk, expected, tolerance = 1e-10)
})
