# ============================================================
# Tests for estimate_values.R: estimate_values()
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# Helper -----------------------------------------------------------------------
make_ev_inputs <- function(augmented = TRUE) {
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

  q_list <- if (augmented) {
    list(
      modelObj::buildModelObj(model = ~ a1, solver.method = "lm",
                              predict.method = "predict.lm"),
      modelObj::buildModelObj(model = ~ a2, solver.method = "lm",
                              predict.method = "predict.lm")
    )
  } else {
    NULL
  }

  list(df = df, nus = nus, pis = pis, regime_all = regime_all,
       q_list = q_list, K = K)
}

# Return structure ------------------------------------------------------------

test_that("estimate_values returns a list with value, df, and q_all", {
  inp <- make_ev_inputs()
  out <- estimate_values(inp$df, inp$q_list, inp$regime_all, TRUE, inp$pis, inp$nus)
  expect_type(out, "list")
  expect_named(out, c("value", "df", "q_all"))
})

# value component -------------------------------------------------------------

test_that("estimate_values value has length equal to number of regimes", {
  inp <- make_ev_inputs()
  out <- estimate_values(inp$df, inp$q_list, inp$regime_all, TRUE, inp$pis, inp$nus)
  expect_length(out$value, length(inp$regime_all))
})

test_that("estimate_values value is numeric and finite", {
  inp <- make_ev_inputs()
  out <- estimate_values(inp$df, inp$q_list, inp$regime_all, TRUE, inp$pis, inp$nus)
  expect_true(is.numeric(out$value))
  expect_true(all(is.finite(out$value)))
})

test_that("estimate_values value equals sum(vTerms) / ns for each regime", {
  inp <- make_ev_inputs()
  out <- estimate_values(inp$df, inp$q_list, inp$regime_all, TRUE, inp$pis, inp$nus)
  ns  <- inp$nus$ns
  for (ell in seq_along(out$value)) {
    expected <- sum(out$df[[ell]]) / ns
    expect_equal(out$value[[ell]], expected, tolerance = 1e-10,
                 label = paste("regime", ell))
  }
})

# df component ----------------------------------------------------------------

test_that("estimate_values df is a list of length L", {
  inp <- make_ev_inputs()
  out <- estimate_values(inp$df, inp$q_list, inp$regime_all, TRUE, inp$pis, inp$nus)
  expect_type(out$df, "list")
  expect_length(out$df, length(inp$regime_all))
})

test_that("estimate_values df matrices have nrow(df) rows and 2K+1 columns", {
  inp <- make_ev_inputs()
  out <- estimate_values(inp$df, inp$q_list, inp$regime_all, TRUE, inp$pis, inp$nus)
  for (ell in seq_along(out$df)) {
    expect_equal(nrow(out$df[[ell]]), nrow(inp$df),
                 label = paste("regime", ell, "row count"))
    expect_equal(ncol(out$df[[ell]]), 2L * inp$K + 1L,
                 label = paste("regime", ell, "column count"))
  }
})

test_that("estimate_values df matrices contain only finite values (augmented)", {
  inp <- make_ev_inputs(augmented = TRUE)
  out <- estimate_values(inp$df, inp$q_list, inp$regime_all, TRUE, inp$pis, inp$nus)
  for (ell in seq_along(out$df)) {
    expect_true(all(is.finite(out$df[[ell]])),
                label = paste("regime", ell))
  }
})

# q_all component -------------------------------------------------------------

test_that("estimate_values q_all is empty when q_list is NULL (IPW)", {
  inp <- make_ev_inputs(augmented = FALSE)
  out <- estimate_values(inp$df, NULL, inp$regime_all, TRUE, inp$pis, inp$nus)
  expect_length(out$q_all, 0L)
})

test_that("estimate_values q_all has length L when q_list is provided (augmented)", {
  inp <- make_ev_inputs(augmented = TRUE)
  out <- estimate_values(inp$df, inp$q_list, inp$regime_all, TRUE, inp$pis, inp$nus)
  expect_length(out$q_all, length(inp$regime_all))
})

# IPW path --------------------------------------------------------------------

test_that("estimate_values IPW values are numeric and finite", {
  inp <- make_ev_inputs(augmented = FALSE)
  out <- estimate_values(inp$df, NULL, inp$regime_all, TRUE, inp$pis, inp$nus)
  expect_true(is.numeric(out$value))
  expect_true(all(is.finite(out$value)))
})

test_that("estimate_values IPW df first 2K columns are all zero", {
  inp <- make_ev_inputs(augmented = FALSE)
  out <- estimate_values(inp$df, NULL, inp$regime_all, TRUE, inp$pis, inp$nus)
  K   <- inp$K
  for (ell in seq_along(out$df)) {
    aug_cols <- out$df[[ell]][, 1:(2L * K), drop = FALSE]
    expect_true(all(aug_cols == 0),
                label = paste("IPW regime", ell, "augmentation cols are zero"))
  }
})
