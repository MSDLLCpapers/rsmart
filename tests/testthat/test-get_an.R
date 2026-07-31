# ============================================================
# Tests for get_an (An.R)
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# Helper -----------------------------------------------------------------------
# Prepares all inputs required by get_an using the pcsttrial dataset.
# augmented = TRUE  -> include Q-function models (AIPWE path)
# augmented = FALSE -> no Q-function models (IPW path)
make_an_inputs <- function(augmented = TRUE) {
  data("pcsttrial", package = "rsmart", envir = environment())

  df    <- pcsttrial
  df$t1 <- df$study_day_enroll
  df$t2 <- df$study_day_rerand
  df$t3 <- df$study_day_outcome
  df$y  <- df$pctchange

  K   <- 2L
  t_s <- max(df$t3)

  df$kappa <- get_kappa(df, t_s, K)
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

  q_list <- if (augmented) {
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
    list(q1, q2)
  } else {
    NULL
  }

  vHats <- estimate_values(df, q_list, regime_all, TRUE, pis, nus)

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

# Tests: helper function ---------------------------------------------------

test_that("make_an_inputs helper returns a well-formed list", {
  inp <- make_an_inputs(augmented = TRUE)
  expect_type(inp, "list")
  expect_named(inp,
               c("df", "nus", "pis", "p_fits", "regime_all",
                 "q_list", "values", "dfs", "q_all"),
               ignore.order = TRUE)
  expect_s3_class(inp$df, "data.frame")
  expect_type(inp$values, "double")
  expect_length(inp$values,     4L)  # 4 embedded regimes
  expect_length(inp$p_fits,     2L)  # K = 2 stages
  expect_length(inp$q_list,     2L)  # K = 2 stages
  expect_length(inp$regime_all, 4L)  # 4 embedded regimes
})

# Tests: IPW path (q_list = NULL) ----------------------------------------------

test_that("get_an returns a matrix for IPW estimation", {
  inp <- make_an_inputs(augmented = FALSE)
  An  <- get_an(inp$df, inp$pis, inp$p_fits, inp$nus, inp$q_all,
               inp$values, inp$regime_all, inp$dfs, TRUE, inp$q_list)
  expect_true(is.matrix(An) || inherits(An, "Matrix"))
})

test_that("get_an has the expected dimension for IPW estimation", {
  inp    <- make_an_inputs(augmented = FALSE)
  An     <- get_an(inp$df, inp$pis, inp$p_fits, inp$nus, inp$q_all,
                  inp$values, inp$regime_all, inp$dfs, TRUE, inp$q_list)
  An_mat <- as.matrix(An)
  # Dimension must equal: (pi params) + (nu params) + (L value params)
  n_pi         <- sum(sapply(inp$p_fits, function(f) length(f@fitObj$coefficients)))
  n_nu         <- length(inp$nus$nu)
  n_L          <- length(inp$values)
  expected_dim <- n_pi + n_nu + n_L
  expect_equal(nrow(An_mat), expected_dim)
  expect_equal(ncol(An_mat), expected_dim)
})

test_that("get_an contains only finite values for IPW estimation", {
  inp <- make_an_inputs(augmented = FALSE)
  An  <- get_an(inp$df, inp$pis, inp$p_fits, inp$nus, inp$q_all,
               inp$values, inp$regime_all, inp$dfs, TRUE, inp$q_list)
  expect_true(all(is.finite(as.vector(An))))
})

test_that("get_an is non-singular (invertible) for IPW estimation", {
  inp    <- make_an_inputs(augmented = FALSE)
  An     <- get_an(inp$df, inp$pis, inp$p_fits, inp$nus, inp$q_all,
                  inp$values, inp$regime_all, inp$dfs, TRUE, inp$q_list)
  An_mat <- as.matrix(An)
  # rcond is the reciprocal condition number; > machine epsilon means invertible
  expect_gt(rcond(An_mat), .Machine$double.eps)
})

test_that("get_an and get_bn have matching dimensions for IPW estimation", {
  inp <- make_an_inputs(augmented = FALSE)
  An  <- get_an(inp$df, inp$pis, inp$p_fits, inp$nus, inp$q_all,
               inp$values, inp$regime_all, inp$dfs, TRUE, inp$q_list)
  Bn  <- get_bn(inp$df, inp$p_fits, inp$nus, inp$q_all,
               inp$dfs, inp$values, TRUE, inp$q_list)
  expect_equal(dim(as.matrix(An)), dim(Bn))
})

# Tests: augmented path (q_list specified) -------------------------------------

test_that("get_an returns a matrix for augmented estimation", {
  inp <- make_an_inputs(augmented = TRUE)
  An  <- get_an(inp$df, inp$pis, inp$p_fits, inp$nus, inp$q_all,
               inp$values, inp$regime_all, inp$dfs, TRUE, inp$q_list)
  expect_true(is.matrix(An) || inherits(An, "Matrix"))
})

test_that("get_an has the expected dimension for augmented estimation", {
  inp    <- make_an_inputs(augmented = TRUE)
  An     <- get_an(inp$df, inp$pis, inp$p_fits, inp$nus, inp$q_all,
                  inp$values, inp$regime_all, inp$dfs, TRUE, inp$q_list)
  An_mat <- as.matrix(An)
  # Dimension must equal: (pi params) + (nu params) + (beta params across all regimes) + (L value params)
  n_pi         <- sum(sapply(inp$p_fits, function(f) length(f@fitObj$coefficients)))
  n_nu         <- length(inp$nus$nu)
  n_beta       <- ncol(ee_psi_beta(inp$df, inp$q_all, TRUE, inp$q_list))
  n_L          <- length(inp$values)
  expected_dim <- n_pi + n_nu + n_beta + n_L
  expect_equal(nrow(An_mat), expected_dim)
  expect_equal(ncol(An_mat), expected_dim)
})

test_that("get_an contains only finite values for augmented estimation", {
  inp <- make_an_inputs(augmented = TRUE)
  An  <- get_an(inp$df, inp$pis, inp$p_fits, inp$nus, inp$q_all,
               inp$values, inp$regime_all, inp$dfs, TRUE, inp$q_list)
  expect_true(all(is.finite(as.vector(An))))
})

test_that("get_an is non-singular (invertible) for augmented estimation", {
  inp    <- make_an_inputs(augmented = TRUE)
  An     <- get_an(inp$df, inp$pis, inp$p_fits, inp$nus, inp$q_all,
                  inp$values, inp$regime_all, inp$dfs, TRUE, inp$q_list)
  An_mat <- as.matrix(An)
  expect_gt(rcond(An_mat), .Machine$double.eps)
})

test_that("get_an and get_bn have matching dimensions for augmented estimation", {
  inp <- make_an_inputs(augmented = TRUE)
  An  <- get_an(inp$df, inp$pis, inp$p_fits, inp$nus, inp$q_all,
               inp$values, inp$regime_all, inp$dfs, TRUE, inp$q_list)
  Bn  <- get_bn(inp$df, inp$p_fits, inp$nus, inp$q_all,
               inp$dfs, inp$values, TRUE, inp$q_list)
  expect_equal(dim(as.matrix(An)), dim(Bn))
})

# Tests: structural differences between IPW and augmented ----------------------

test_that("get_an augmented matrix is larger than IPW matrix", {
  inp_ipw <- make_an_inputs(augmented = FALSE)
  inp_aug <- make_an_inputs(augmented = TRUE)

  An_ipw  <- get_an(inp_ipw$df, inp_ipw$pis, inp_ipw$p_fits, inp_ipw$nus,
                   inp_ipw$q_all, inp_ipw$values, inp_ipw$regime_all,
                   inp_ipw$dfs, TRUE, inp_ipw$q_list)
  An_aug  <- get_an(inp_aug$df, inp_aug$pis, inp_aug$p_fits, inp_aug$nus,
                   inp_aug$q_all, inp_aug$values, inp_aug$regime_all,
                   inp_aug$dfs, TRUE, inp_aug$q_list)

  expect_gt(nrow(as.matrix(An_aug)), nrow(as.matrix(An_ipw)))
})
