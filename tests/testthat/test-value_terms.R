# ============================================================
# Tests for value_terms.R: value_terms()
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# Helper -----------------------------------------------------------------------
make_vterms_inputs <- function(augmented = TRUE) {
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
  pis <- pi_fits(df, list(p1, p2))[["ps"]]

  regime_all <- regime_list_no_trt_resp(
    emb_regimes = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
    dat        = df,
    resp_trt    = list("r2" = list(0, 0, 0, 0))
  )

  # Use regime 1 for single-regime tests
  regime1    <- regime_all[[1]]
  regime_ind <- regime1$regime_ind
  regime_mat <- regime1$regime

  qs <- if (augmented) {
    q1 <- modelObj::buildModelObj(model = ~ a1, solver.method = "lm",
                                  predict.method = "predict.lm")
    q2 <- modelObj::buildModelObj(model = ~ a2, solver.method = "lm",
                                  predict.method = "predict.lm")
    get_q_fits(df, list(q1, q2), regime_mat, feasible_sets_indicator = TRUE)
  } else {
    NULL
  }

  list(df = df, nus = nus, pis = pis, regime_ind = regime_ind, qs = qs, K = K)
}

# Return dimensions ------------------------------------------------------------

test_that("value_terms returns a matrix with nrow(df) rows and 2K+1 columns (augmented)", {
  inp <- make_vterms_inputs(augmented = TRUE)
  out <- value_terms(inp$df, inp$regime_ind, inp$pis, inp$qs, inp$nus)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), nrow(inp$df))
  expect_equal(ncol(out), 2L * inp$K + 1L)
})

test_that("value_terms returns a matrix with nrow(df) rows and 2K+1 columns (IPW)", {
  inp <- make_vterms_inputs(augmented = FALSE)
  out <- value_terms(inp$df, inp$regime_ind, inp$pis, inp$qs, inp$nus)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), nrow(inp$df))
  expect_equal(ncol(out), 2L * inp$K + 1L)
})

# Finite values ----------------------------------------------------------------

test_that("value_terms contains only finite values (augmented)", {
  inp <- make_vterms_inputs(augmented = TRUE)
  out <- value_terms(inp$df, inp$regime_ind, inp$pis, inp$qs, inp$nus)
  expect_true(all(is.finite(out)))
})

test_that("value_terms contains only finite values (IPW)", {
  inp <- make_vterms_inputs(augmented = FALSE)
  out <- value_terms(inp$df, inp$regime_ind, inp$pis, inp$qs, inp$nus)
  expect_true(all(is.finite(out)))
})

# IPW-specific behaviour -------------------------------------------------------

test_that("value_terms IPW: first 2K augmentation columns are all zero", {
  inp <- make_vterms_inputs(augmented = FALSE)
  out <- value_terms(inp$df, inp$regime_ind, inp$pis, inp$qs, inp$nus)
  aug_cols <- out[, 1:(2L * inp$K), drop = FALSE]
  expect_true(all(aug_cols == 0))
})

test_that("value_terms IPW: last column is non-zero for subjects with R = 2K+1", {
  inp <- make_vterms_inputs(augmented = FALSE)
  out <- value_terms(inp$df, inp$regime_ind, inp$pis, inp$qs, inp$nus)
  # At least some individuals should have reached the final coarsening level
  expect_true(any(out[, 2L * inp$K + 1L] != 0))
})

# Augmented-specific behaviour -------------------------------------------------

test_that("value_terms augmented: some augmentation columns are non-zero", {
  inp <- make_vterms_inputs(augmented = TRUE)
  out <- value_terms(inp$df, inp$regime_ind, inp$pis, inp$qs, inp$nus)
  aug_cols <- out[, 1:(2L * inp$K), drop = FALSE]
  expect_true(any(aug_cols != 0))
})

# Value consistency ------------------------------------------------------------

test_that("value_terms sum / ns equals the estimate_values regime value", {
  inp <- make_vterms_inputs(augmented = TRUE)
  out <- value_terms(inp$df, inp$regime_ind, inp$pis, inp$qs, inp$nus)
  # Compute the value the same way estimate_values does
  value_from_vterms <- sum(out) / inp$nus$ns
  # Cross-check against estimate_values for regime 1
  regime_all <- regime_list_no_trt_resp(
    emb_regimes = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
    dat        = inp$df,
    resp_trt    = list("r2" = list(0, 0, 0, 0))
  )
  q1 <- modelObj::buildModelObj(model = ~ a1, solver.method = "lm",
                                predict.method = "predict.lm")
  q2 <- modelObj::buildModelObj(model = ~ a2, solver.method = "lm",
                                predict.method = "predict.lm")
  ev <- estimate_values(inp$df, list(q1, q2), regime_all, TRUE, inp$pis, inp$nus)
  expect_equal(value_from_vterms, ev$value[[1]], tolerance = 1e-10)
})

test_that("value_terms IPW sum / ns equals the estimate_values IPW regime value", {
  inp <- make_vterms_inputs(augmented = FALSE)
  out <- value_terms(inp$df, inp$regime_ind, inp$pis, inp$qs, inp$nus)
  value_from_vterms <- sum(out) / inp$nus$ns

  regime_all <- regime_list_no_trt_resp(
    emb_regimes = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
    dat        = inp$df,
    resp_trt    = list("r2" = list(0, 0, 0, 0))
  )
  ev <- estimate_values(inp$df, NULL, regime_all, TRUE, inp$pis, inp$nus)
  expect_equal(value_from_vterms, ev$value[[1]], tolerance = 1e-10)
})
