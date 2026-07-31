# ============================================================
# Tests for block_rand
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

test_that("block_rand returns a function", {
  fn <- block_rand(block_rep = 2)
  expect_true(is.function(fn))
})

test_that("block_rand returned function has correct signature", {
  fn <- block_rand(block_rep = 2)
  args <- names(formals(fn))
  expect_equal(args, c("n", "n_treatments", "prob"))
})

test_that("block_rand returns correct length vector", {
  fn <- block_rand(block_rep = 2)
  set.seed(1)
  result <- fn(n = 50, n_treatments = 2, prob = c(0.5, 0.5))
  expect_length(result, 50)
})

test_that("block_rand returns valid treatment codes for 2 treatments", {
  fn <- block_rand(block_rep = 2)
  set.seed(1)
  result <- fn(n = 100, n_treatments = 2, prob = c(0.5, 0.5))
  expect_true(all(result %in% c(0, 1)))
})

test_that("block_rand returns valid treatment codes for 3 treatments", {
  fn <- block_rand(block_rep = 2)
  set.seed(1)
  result <- fn(n = 300, n_treatments = 3, prob = c(1/3, 1/3, 1/3))
  expect_true(all(result %in% c(0, 1, 2)))
})

test_that("block_rand gives exactly equal counts when n is a multiple of block size", {
  fn <- block_rand(block_rep = 2)
  # block size = 2 * 2 = 4; n = 100 is 25 blocks exactly
  set.seed(1)
  result <- fn(n = 100, n_treatments = 2, prob = c(0.5, 0.5))
  counts <- table(result)
  expect_equal(as.integer(counts["0"]), 50)
  expect_equal(as.integer(counts["1"]), 50)
})

test_that("block_rand gives exactly equal counts for 3 treatments when n is a multiple", {
  fn <- block_rand(block_rep = 3)
  # block size = 3 * 3 = 9; n = 90 is 10 blocks exactly
  set.seed(1)
  result <- fn(n = 90, n_treatments = 3, prob = c(1/3, 1/3, 1/3))
  counts <- table(result)
  expect_equal(as.integer(counts["0"]), 30)
  expect_equal(as.integer(counts["1"]), 30)
  expect_equal(as.integer(counts["2"]), 30)
})

test_that("block_rand gives near-equal counts when n is not a multiple of block size", {
  fn <- block_rand(block_rep = 2)
  # block size = 4; n = 101 -> 26 blocks generated, truncated to 101
  set.seed(1)
  result <- fn(n = 101, n_treatments = 2, prob = c(0.5, 0.5))
  counts <- table(result)
  # difference should be at most 1 (one incomplete block)
  expect_true(abs(counts["0"] - counts["1"]) <= 1)
})

test_that("block_rand has balance within each complete block", {
  fn <- block_rand(block_rep = 2)
  set.seed(42)
  # 5 complete blocks of size 4
  result <- fn(n = 20, n_treatments = 2, prob = c(0.5, 0.5))
  for (b in 1:5) {
    block <- result[((b - 1) * 4 + 1):(b * 4)]
    expect_equal(sum(block == 0), 2)
    expect_equal(sum(block == 1), 2)
  }
})

test_that("block_rand with block_rep = 1 produces blocks of size n_treatments", {
  fn <- block_rand(block_rep = 1)
  set.seed(1)
  result <- fn(n = 20, n_treatments = 2, prob = c(0.5, 0.5))
  # 10 blocks of size 2, each with one 0 and one 1
  for (b in 1:10) {
    block <- result[((b - 1) * 2 + 1):(b * 2)]
    expect_equal(sort(block), c(0, 1))
  }
})

test_that("block_rand ignores the prob argument", {
  fn <- block_rand(block_rep = 2)
  set.seed(1)
  result1 <- fn(n = 100, n_treatments = 2, prob = c(0.5, 0.5))
  set.seed(1)
  result2 <- fn(n = 100, n_treatments = 2, prob = c(0.9, 0.1))
  expect_identical(result1, result2)
})

test_that("block_rand is deterministic with the same seed", {
  fn <- block_rand(block_rep = 2)
  set.seed(123)
  result1 <- fn(n = 50, n_treatments = 2, prob = c(0.5, 0.5))
  set.seed(123)
  result2 <- fn(n = 50, n_treatments = 2, prob = c(0.5, 0.5))
  expect_identical(result1, result2)
})

test_that("block_rand stops for invalid block_rep", {
  expect_error(block_rand(block_rep = 0), "block_rep must be a positive integer")
  expect_error(block_rand(block_rep = -1), "block_rep must be a positive integer")
  expect_error(block_rand(block_rep = 1.5), "block_rep must be a positive integer")
  expect_error(block_rand(block_rep = "a"), "block_rep must be a positive integer")
})

test_that("block_rand works with sim_treatment", {
  set.seed(1)
  result <- sim_treatment(n = 100, n_treatments = 2,
                          rand_prob_fn = block_rand(block_rep = 2))
  expect_length(result, 100)
  expect_true(all(result %in% c(0L, 1L)))
  # should be exactly balanced
  expect_equal(sum(result == 0), 50)
  expect_equal(sum(result == 1), 50)
})
