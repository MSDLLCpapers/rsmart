# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))


# Helper to set up a small SMART dataset and model inputs
setup_iaipwe_inputs <- function(n = 200, seed = 42) {
  set.seed(seed)
  df <- gen_no_trt_resp(n = n, s2 = 100, block_rep = 2, r2p = 0.5)

  p1 <- modelObj::buildModelObj(
    model = ~1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  p2 <- modelObj::buildModelObj(
    model = ~ I(a1 == 0):I(r2 == 0) - 1,
    solver.method = "glm",
    solver.args = list(family = "binomial"),
    predict.method = "predict.glm",
    predict.args = list(type = "response")
  )

  pi_list <- list(p1, p2)

  q1 <- modelObj::buildModelObj(
    model = ~ x11 + x12 + a1 + x11:a1,
    solver.method = "lm",
    predict.method = "predict.lm"
  )

  q2 <- modelObj::buildModelObj(
    model = ~ x11 + x12 + x21 + a1 + a2 + a1:a2,
    solver.method = "lm",
    predict.method = "predict.lm"
  )

  q_list <- list(q1, q2)

  regime_all <- regime_list_no_trt_resp(
    emb_regimes = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
    df,
    resp_trt = list("r2" = list(0, 0, 0, 0))
  )

  t_s <- max(df$t3)

  list(
    df = df, pi_list = pi_list, q_list = q_list,
    regime_all = regime_all, feasible_sets_indicator = TRUE, t_s = t_s
  )
}

test_that("asymptotic and bootstrap return lists have the same names", {
  inputs <- setup_iaipwe_inputs()

  res_asym <- iaipwe(df = inputs$df,
                     pi_list = inputs$pi_list,
                     q_list = inputs$q_list,
                     regime_all = inputs$regime_all,
                     feasible_sets_indicator = TRUE,
                     t_s = inputs$t_s,
                     B = NULL)

  res_boot <- iaipwe(df = inputs$df,
                     pi_list = inputs$pi_list,
                     q_list = inputs$q_list,
                     regime_all = inputs$regime_all,
                     feasible_sets_indicator = TRUE,
                     t_s = inputs$t_s,
                     B = 50)

  expect_equal(names(res_asym), names(res_boot))
})

test_that("covariance is L x L for both asymptotic and bootstrap", {
  inputs <- setup_iaipwe_inputs()
  L <- length(inputs$regime_all)

  res_asym <- iaipwe(df = inputs$df,
                     pi_list = inputs$pi_list,
                     q_list = inputs$q_list,
                     regime_all = inputs$regime_all,
                     feasible_sets_indicator = TRUE,
                     t_s = inputs$t_s,
                     B = NULL)

  res_boot <- iaipwe(df = inputs$df,
                     pi_list = inputs$pi_list,
                     q_list = inputs$q_list,
                     regime_all = inputs$regime_all,
                     feasible_sets_indicator = TRUE,
                     t_s = inputs$t_s,
                     B = 50)

  expect_equal(dim(res_asym$covariance), c(L, L))
  expect_equal(dim(res_boot$covariance), c(L, L))
})

test_that("values and se have the same length for both methods", {
  inputs <- setup_iaipwe_inputs()
  L <- length(inputs$regime_all)

  res_asym <- iaipwe(df = inputs$df,
                     pi_list = inputs$pi_list,
                     q_list = inputs$q_list,
                     regime_all = inputs$regime_all,
                     feasible_sets_indicator = TRUE,
                     t_s = inputs$t_s,
                     B = NULL)

  res_boot <- iaipwe(df = inputs$df,
                     pi_list = inputs$pi_list,
                     q_list = inputs$q_list,
                     regime_all = inputs$regime_all,
                     feasible_sets_indicator = TRUE,
                     t_s = inputs$t_s,
                     B = 50)

  # values should be the same point estimates regardless of variance method
  expect_equal(res_asym$values, res_boot$values)

  # se should have length L in both cases
  expect_length(res_asym$se, L)
  expect_length(res_boot$se, L)
})

test_that("variance_choice is correctly labeled", {
  inputs <- setup_iaipwe_inputs()

  res_asym <- iaipwe(df = inputs$df,
                     pi_list = inputs$pi_list,
                     q_list = inputs$q_list,
                     regime_all = inputs$regime_all,
                     feasible_sets_indicator = TRUE,
                     t_s = inputs$t_s,
                     B = NULL)

  res_boot <- iaipwe(df = inputs$df,
                     pi_list = inputs$pi_list,
                     q_list = inputs$q_list,
                     regime_all = inputs$regime_all,
                     feasible_sets_indicator = TRUE,
                     t_s = inputs$t_s,
                     B = 50)

  expect_equal(res_asym$variance_choice, "Asymptotic")
  expect_equal(res_boot$variance_choice, "Bootstrap")
})

test_that("chi_square structure is consistent between methods", {
  inputs <- setup_iaipwe_inputs()

  res_asym <- iaipwe(df = inputs$df,
                     pi_list = inputs$pi_list,
                     q_list = inputs$q_list,
                     regime_all = inputs$regime_all,
                     feasible_sets_indicator = TRUE,
                     t_s = inputs$t_s,
                     B = NULL)

  res_boot <- iaipwe(df = inputs$df,
                     pi_list = inputs$pi_list,
                     q_list = inputs$q_list,
                     regime_all = inputs$regime_all,
                     feasible_sets_indicator = TRUE,
                     t_s = inputs$t_s,
                     B = 50)

  expect_named(res_asym$chi_square, c("Statistic", "p_value"))
  expect_named(res_boot$chi_square, c("Statistic", "p_value"))

  # Both should be numeric scalars
  expect_length(res_asym$chi_square$Statistic, 1)
  expect_length(res_boot$chi_square$Statistic, 1)
  expect_true(res_asym$chi_square$p_value >= 0 && res_asym$chi_square$p_value <= 1)
  expect_true(res_boot$chi_square$p_value >= 0 && res_boot$chi_square$p_value <= 1)
})

test_that("bootstrap and asymptotic covariances are on the same scale", {
  inputs <- setup_iaipwe_inputs()
  set.seed(123)
  result_asymp <- do.call(iaipwe, c(inputs, list(B = NULL)))
  set.seed(123)
  result_boot <- do.call(iaipwe, c(inputs, list(B = 50)))

  # Both covariances should be of similar magnitude
  # (not off by a factor of n)
  ratio <- mean(diag(result_boot$covariance)) / mean(diag(result_asymp$covariance))
  expect_equal(ratio, 1, tolerance = 0.5)  # Allow some tolerance due to randomness
})
