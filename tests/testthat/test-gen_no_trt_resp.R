# ============================================================
# Tests for gen_no_trt_resp.R:
#   gen_no_trt_resp() and regime_list_no_trt_resp()
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# ---- gen_no_trt_resp --------------------------------------------------------

test_that("gen_no_trt_resp returns a data frame with n rows", {
  set.seed(1)
  dat <- gen_no_trt_resp(n = 100)
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 100L)
})

test_that("gen_no_trt_resp returns required columns", {
  set.seed(1)
  dat <- gen_no_trt_resp(n = 100)
  required <- c("t1", "t2", "t3", "a1", "r2", "a2", "x11", "x12", "x21", "y", "id")
  expect_true(all(required %in% colnames(dat)))
})

test_that("gen_no_trt_resp stage-1 treatment a1 is binary", {
  set.seed(1)
  dat <- gen_no_trt_resp(n = 200)
  expect_true(all(dat$a1 %in% c(0L, 1L)))
})

test_that("gen_no_trt_resp stage-2 treatment a2 is binary", {
  set.seed(1)
  dat <- gen_no_trt_resp(n = 200)
  expect_true(all(dat$a2 %in% c(0L, 1L)))
})

test_that("gen_no_trt_resp response indicator r2 is binary", {
  set.seed(1)
  dat <- gen_no_trt_resp(n = 200)
  expect_true(all(dat$r2 %in% c(0L, 1L)))
})

test_that("gen_no_trt_resp responders have a2 = 0 (not re-randomized)", {
  set.seed(1)
  dat <- gen_no_trt_resp(n = 400, r2p = 0.5)
  responders <- dat[dat$r2 == 1, ]
  expect_true(all(responders$a2 == 0L))
})

test_that("gen_no_trt_resp timing satisfies t1 <= t2 <= t3", {
  set.seed(1)
  dat <- gen_no_trt_resp(n = 200)
  expect_true(all(dat$t1 <= dat$t2))
  expect_true(all(dat$t2 <= dat$t3))
})

test_that("gen_no_trt_resp id equals 1:n", {
  set.seed(1)
  dat <- gen_no_trt_resp(n = 50)
  expect_equal(sort(dat$id), 1:50)
})

test_that("gen_no_trt_resp y is numeric and finite", {
  set.seed(1)
  dat <- gen_no_trt_resp(n = 100)
  expect_true(is.numeric(dat$y))
  expect_true(all(is.finite(dat$y)))
})

test_that("gen_no_trt_resp is reproducible with set.seed", {
  set.seed(42)
  dat1 <- gen_no_trt_resp(n = 50)
  set.seed(42)
  dat2 <- gen_no_trt_resp(n = 50)
  expect_equal(dat1, dat2)
})

# Error handling ---------------------------------------------------------------

test_that("gen_no_trt_resp errors on negative n", {
  expect_error(gen_no_trt_resp(n = -1), regexp = "non-negative integer")
})

test_that("gen_no_trt_resp errors on non-integer n", {
  expect_error(gen_no_trt_resp(n = 10.5), regexp = "non-negative integer")
})

test_that("gen_no_trt_resp errors on non-positive s2", {
  expect_error(gen_no_trt_resp(n = 100, s2 = 0), regexp = "positive numeric")
  expect_error(gen_no_trt_resp(n = 100, s2 = -5), regexp = "positive numeric")
})

test_that("gen_no_trt_resp errors when r2p is outside [0, 1]", {
  expect_error(gen_no_trt_resp(n = 100, r2p = -0.1), regexp = "between 0 and 1")
  expect_error(gen_no_trt_resp(n = 100, r2p = 1.1),  regexp = "between 0 and 1")
})

test_that("gen_no_trt_resp errors on non-integer block_rep", {
  expect_error(gen_no_trt_resp(n = 100, block_rep = 1.5), regexp = "non-negative integer")
})

# ---- regime_list_no_trt_resp -------------------------------------------------

test_that("regime_list_no_trt_resp returns a list with length equal to emb_regimes", {
  set.seed(1)
  dat  <- gen_no_trt_resp(n = 100)
  regs <- regime_list_no_trt_resp(
    emb_regimes = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
    dat        = dat,
    resp_trt    = list("r2" = list(0, 0, 0, 0))
  )
  expect_type(regs, "list")
  expect_length(regs, 4L)
})

test_that("regime_list_no_trt_resp each element has regime and regime_ind", {
  set.seed(1)
  dat  <- gen_no_trt_resp(n = 100)
  regs <- regime_list_no_trt_resp(
    emb_regimes = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
    dat        = dat,
    resp_trt    = list("r2" = list(0, 0, 0, 0))
  )
  for (ell in seq_along(regs)) {
    expect_named(regs[[ell]], c("regime", "regime_ind"),
                 label = paste("regime", ell))
  }
})

test_that("regime_list_no_trt_resp regime matrix has n rows and K columns", {
  set.seed(1)
  n   <- 100L
  K   <- 2L
  dat <- gen_no_trt_resp(n = n)
  regs <- regime_list_no_trt_resp(
    emb_regimes = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
    dat        = dat,
    resp_trt    = list("r2" = list(0, 0, 0, 0))
  )
  for (ell in seq_along(regs)) {
    expect_equal(dim(regs[[ell]]$regime),     c(n, K), label = paste("regime",     ell))
    expect_equal(dim(regs[[ell]]$regime_ind), c(n, K), label = paste("regime_ind", ell))
  }
})

test_that("regime_list_no_trt_resp regime_ind values are in {0, 1}", {
  set.seed(1)
  dat  <- gen_no_trt_resp(n = 200)
  regs <- regime_list_no_trt_resp(
    emb_regimes = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
    dat        = dat,
    resp_trt    = list("r2" = list(0, 0, 0, 0))
  )
  for (ell in seq_along(regs)) {
    expect_true(all(regs[[ell]]$regime_ind %in% c(0L, 1L)),
                label = paste("regime_ind", ell))
  }
})

test_that("regime_list_no_trt_resp responders get the designated responder treatment", {
  set.seed(1)
  dat  <- gen_no_trt_resp(n = 400, r2p = 0.5)
  # All regimes assign responder treatment 0
  regs <- regime_list_no_trt_resp(
    emb_regimes = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1)),
    dat        = dat,
    resp_trt    = list("r2" = list(0, 0, 0, 0))
  )
  responders <- dat$r2 == 1
  for (ell in seq_along(regs)) {
    expect_true(all(regs[[ell]]$regime[responders, 2] == 0),
                label = paste("responder treatment in regime", ell))
  }
})
