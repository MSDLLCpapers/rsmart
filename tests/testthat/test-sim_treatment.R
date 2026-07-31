# ============================================================
# Tests for sim_treatment
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# --- Return type and structure ---

test_that("sim_treatment returns integer vector when no dat provided", {
  set.seed(1)
  result <- sim_treatment(n = 50, n_treatments = 2)
  expect_type(result, "integer")
  expect_length(result, 50)
})

test_that("sim_treatment returns data frame with new column when dat provided", {
  set.seed(1)
  df <- data.frame(x1 = rnorm(50))
  result <- sim_treatment(n_treatments = 2, dat = df)
  expect_s3_class(result, "data.frame")
  expect_true("a1" %in% colnames(result))
  expect_equal(nrow(result), 50)
})

test_that("sim_treatment column is named a<stage>", {
  set.seed(1)
  df <- data.frame(a1 = rbinom(50, 1, 0.5))
  result <- sim_treatment(n_treatments = 2, dat = df, stage = 2)
  expect_true("a2" %in% colnames(result))
})

# --- Treatment codes ---

test_that("sim_treatment returns valid treatment codes for 2 treatments", {
  set.seed(1)
  result <- sim_treatment(n = 200, n_treatments = 2)
  expect_true(all(result %in% c(0L, 1L)))
})

test_that("sim_treatment returns valid treatment codes for 3 treatments", {
  set.seed(1)
  result <- sim_treatment(n = 300, n_treatments = 3)
  expect_true(all(result %in% c(0L, 1L, 2L)))
})

# --- prob argument ---

test_that("sim_treatment with default prob gives approximately equal split", {
  set.seed(42)
  result <- sim_treatment(n = 10000, n_treatments = 2)
  prop <- mean(result == 1)
  expect_equal(prop, 0.5, tolerance = 0.03)
})

test_that("sim_treatment with unequal prob vector shifts allocation", {
  set.seed(42)
  result <- sim_treatment(n = 10000, n_treatments = 2, prob = c(0.8, 0.2))
  prop <- mean(result == 0)
  expect_equal(prop, 0.8, tolerance = 0.03)
})

test_that("sim_treatment stops when prob vector has wrong length", {
  expect_error(
    sim_treatment(n = 50, n_treatments = 2, prob = c(0.5, 0.3, 0.2)),
    "prob vector must have length equal to n_treatments"
  )
})

test_that("sim_treatment stops when prob is not numeric, list, or NULL", {
  expect_error(
    sim_treatment(n = 50, n_treatments = 2, prob = "bad"),
    "prob must be NULL, a numeric vector, or a named list"
  )
})

test_that("sim_treatment with prob as named list uses group-specific probabilities", {
  set.seed(42)
  n <- 10000
  df <- data.frame(a1 = c(rep(0, n / 2), rep(1, n / 2)))
  result <- sim_treatment(
    n_treatments = 2, dat = df, stage = 2,
    prob = list("0" = c(0.9, 0.1), "1" = c(0.1, 0.9))
  )
  # a1==0 group should be ~90% treatment 0
  prop_a1_0 <- mean(result$a2[result$a1 == 0] == 0)
  expect_equal(prop_a1_0, 0.9, tolerance = 0.03)
  # a1==1 group should be ~90% treatment 1
  prop_a1_1 <- mean(result$a2[result$a1 == 1] == 1)
  expect_equal(prop_a1_1, 0.9, tolerance = 0.03)
})

test_that("sim_treatment stops when prob list is missing a group", {
  df <- data.frame(a1 = c(0, 0, 1, 1))
  expect_error(
    sim_treatment(n_treatments = 2, dat = df, stage = 2,
                  prob = list("0" = c(0.5, 0.5))),
    "prob list is missing an entry for group"
  )
})

# --- Input validation ---

test_that("sim_treatment stops for invalid n", {
  expect_error(sim_treatment(n = -1, n_treatments = 2), "n must be a positive integer")
  expect_error(sim_treatment(n = 0, n_treatments = 2), "n must be a positive integer")
  expect_error(sim_treatment(n = 1.5, n_treatments = 2), "n must be a positive integer")
  expect_error(sim_treatment(n_treatments = 2), "n must be a positive integer")
})

test_that("sim_treatment stops for invalid n_treatments", {
  expect_error(sim_treatment(n = 10, n_treatments = 1), "n_treatments must be an integer >= 2")
  expect_error(sim_treatment(n = 10, n_treatments = 2.5), "n_treatments must be an integer >= 2")
})

test_that("sim_treatment stops for invalid stage", {
  expect_error(
    sim_treatment(n = 10, n_treatments = 2, stage = 0),
    "stage must be a positive integer"
  )
  expect_error(
    sim_treatment(n = 10, n_treatments = 2, stage = 1.5),
    "stage must be a positive integer"
  )
})

test_that("sim_treatment stops when dat already contains the target column", {
  df <- data.frame(a1 = c(0, 1))
  expect_error(
    sim_treatment(n_treatments = 2, dat = df, stage = 1),
    "dat already contains a column named 'a1'"
  )
})

test_that("sim_treatment stops when dat is not a data frame", {
  expect_error(
    sim_treatment(n_treatments = 2, dat = list(a1 = c(0, 1))),
    "dat must be a data frame"
  )
})

# --- stage > 1 ---

test_that("sim_treatment stops when stage > 1 and dat is NULL", {
  expect_error(
    sim_treatment(n = 10, n_treatments = 2, stage = 2),
    "dat must be provided when stage > 1"
  )
})

test_that("sim_treatment stops when prior treatment columns are missing", {
  df <- data.frame(x1 = rnorm(10))
  expect_error(
    sim_treatment(n_treatments = 2, dat = df, stage = 2),
    "dat is missing prior treatment columns: a1"
  )
})

test_that("sim_treatment at stage 2 randomizes within a1 groups", {
  # With balanced data and default prob, each a1 group should be ~50/50
  set.seed(42)
  n <- 10000
  df <- data.frame(a1 = c(rep(0, n / 2), rep(1, n / 2)))
  result <- sim_treatment(n_treatments = 2, dat = df, stage = 2)
  for (a1_val in c(0, 1)) {
    prop <- mean(result$a2[result$a1 == a1_val] == 1)
    expect_equal(prop, 0.5, tolerance = 0.05)
  }
})

test_that("sim_treatment at stage 3 requires a1 and a2", {
  df <- data.frame(a1 = c(0, 1, 0, 1))
  expect_error(
    sim_treatment(n_treatments = 2, dat = df, stage = 3),
    "dat is missing prior treatment columns: a2"
  )
})

test_that("sim_treatment at stage 3 randomizes within a1 x a2 groups", {
  set.seed(42)
  n_per <- 5000
  df <- data.frame(
    a1 = rep(c(0, 0, 1, 1), each = n_per),
    a2 = rep(c(0, 1, 0, 1), each = n_per)
  )
  result <- sim_treatment(n_treatments = 2, dat = df, stage = 3)
  expect_true("a3" %in% colnames(result))
  # each group should be ~50/50
  for (a1_val in c(0, 1)) {
    for (a2_val in c(0, 1)) {
      sub <- result$a3[result$a1 == a1_val & result$a2 == a2_val]
      expect_equal(mean(sub == 1), 0.5, tolerance = 0.06)
    }
  }
})

# --- randomize_response ---

test_that("sim_treatment stops for invalid randomize_response", {
  expect_error(
    sim_treatment(n = 10, n_treatments = 2, randomize_response = "X"),
    "randomize_response must be NULL, 'Y', or 'N'"
  )
})

test_that("sim_treatment stops when randomize_response used at stage 1", {
  expect_error(
    sim_treatment(n = 10, n_treatments = 2, randomize_response = "Y"),
    "requires stage > 1"
  )
  expect_error(
    sim_treatment(n = 10, n_treatments = 2, randomize_response = "N"),
    "requires stage > 1"
  )
})

test_that("sim_treatment stops when response column is missing for randomize_response Y", {
  df <- data.frame(a1 = c(0, 1, 0, 1))
  expect_error(
    sim_treatment(n_treatments = 2, dat = df, stage = 2, randomize_response = "Y"),
    "dat is missing response column 'r2'"
  )
})

test_that("sim_treatment stops when response column is missing for randomize_response N", {
  df <- data.frame(a1 = c(0, 1, 0, 1))
  expect_error(
    sim_treatment(n_treatments = 2, dat = df, stage = 2, randomize_response = "N"),
    "dat is missing response column 'r2'"
  )
})

test_that("randomize_response N assigns treatment 0 to all responders", {
  set.seed(42)
  n <- 200
  df <- data.frame(
    a1 = rbinom(n, 1, 0.5),
    r2 = rbinom(n, 1, 0.5)
  )
  result <- sim_treatment(n_treatments = 2, dat = df, stage = 2,
                          randomize_response = "N")
  # all responders should have a2 == 0

  expect_true(all(result$a2[result$r2 == 1] == 0))
})

test_that("randomize_response N randomizes non-responders", {
  set.seed(42)
  n <- 5000
  df <- data.frame(
    a1 = rbinom(n, 1, 0.5),
    r2 = rbinom(n, 1, 0.5)
  )
  result <- sim_treatment(n_treatments = 2, dat = df, stage = 2,
                          randomize_response = "N")
  # non-responders should be approximately 50/50
  non_resp <- result$a2[result$r2 == 0]
  expect_equal(mean(non_resp == 1), 0.5, tolerance = 0.05)
})

test_that("randomize_response Y randomizes within response groups", {
  set.seed(42)
  n <- 5000
  df <- data.frame(
    a1 = rbinom(n, 1, 0.5),
    r2 = rbinom(n, 1, 0.5)
  )
  result <- sim_treatment(n_treatments = 2, dat = df, stage = 2,
                          randomize_response = "Y")
  # both responders and non-responders should have treatments
  expect_true(any(result$a2[result$r2 == 1] == 1))
  expect_true(any(result$a2[result$r2 == 0] == 1))
  # each subgroup should be ~50/50
  expect_equal(mean(result$a2[result$r2 == 0] == 1), 0.5, tolerance = 0.05)
  expect_equal(mean(result$a2[result$r2 == 1] == 1), 0.5, tolerance = 0.05)
})

# --- custom rand_prob_fn ---

test_that("sim_treatment works with a custom rand_prob_fn", {
  always_zero <- function(n, n_treatments, prob) rep(0L, n)
  set.seed(1)
  result <- sim_treatment(n = 50, n_treatments = 2, rand_prob_fn = always_zero)
  expect_true(all(result == 0L))
})

test_that("sim_treatment with block_rand at stage 2 balances within groups", {
  n_per <- 50
  df <- data.frame(a1 = c(rep(0, n_per), rep(1, n_per)))
  set.seed(1)
  result <- sim_treatment(n_treatments = 2, dat = df, stage = 2,
                          rand_prob_fn = block_rand(block_rep = 2))
  # within each a1 group, block rand should give equal counts
  for (a1_val in c(0, 1)) {
    sub <- result$a2[result$a1 == a1_val]
    expect_equal(sum(sub == 0), sum(sub == 1))
  }
})

# --- Determinism ---

test_that("sim_treatment is deterministic with the same seed", {
  set.seed(99)
  r1 <- sim_treatment(n = 100, n_treatments = 2)
  set.seed(99)
  r2 <- sim_treatment(n = 100, n_treatments = 2)
  expect_identical(r1, r2)
})

test_that("sim_treatment with dat is deterministic with the same seed", {
  df <- data.frame(a1 = c(rep(0, 50), rep(1, 50)))
  set.seed(99)
  r1 <- sim_treatment(n_treatments = 2, dat = df, stage = 2)
  set.seed(99)
  r2 <- sim_treatment(n_treatments = 2, dat = df, stage = 2)
  expect_identical(r1, r2)
})
