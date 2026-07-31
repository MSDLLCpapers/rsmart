# ============================================================
# Tests for ee_nu.R: ee_psi_nu and ee_dpsi_nu
# ============================================================
# Skip this entire file if the operating system is NOT windows
testthat::skip_on_os(c("mac", "linux", "solaris"))

# Helper -----------------------------------------------------------------------
make_nu_inputs <- function() {
  data("pcsttrial", package = "rsmart", envir = environment())
  df       <- pcsttrial
  df$t1    <- df$study_day_enroll
  df$t2    <- df$study_day_rerand
  df$t3    <- df$study_day_outcome
  K        <- 2L
  df$kappa <- get_kappa(df, max(df$t3), K)
  nus      <- get_nu(df, K)
  list(df = df, nus = nus, K = K)
}

# ee_psi_nu --------------------------------------------------------------------

test_that("ee_psi_nu returns a matrix with n rows and K+1 columns", {
  inp <- make_nu_inputs()
  out <- ee_psi_nu(inp$df, inp$nus)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), nrow(inp$df))
  expect_equal(ncol(out), inp$K + 1L)
})

test_that("ee_psi_nu column sums are exactly zero at estimated nu", {
  inp  <- make_nu_inputs()
  out  <- ee_psi_nu(inp$df, inp$nus)
  sums <- colSums(out)
  expect_equal(sums, rep(0, inp$K + 1L), tolerance = 1e-10)
})

test_that("ee_psi_nu contains only finite values", {
  inp <- make_nu_inputs()
  out <- ee_psi_nu(inp$df, inp$nus)
  expect_true(all(is.finite(out)))
})

test_that("ee_psi_nu column k equals I(kappa>=k) minus nu_k * I(kappa>=1)", {
  inp <- make_nu_inputs()
  out <- ee_psi_nu(inp$df, inp$nus)
  df  <- inp$df
  # spot-check stage 1 (first column)
  expected_col1 <- (df$kappa >= 1) - (df$kappa >= 1) * inp$nus$nu[[1]]
  expect_equal(as.numeric(out[, 1]), as.numeric(expected_col1), tolerance = 1e-12)
})

# ee_dpsi_nu -------------------------------------------------------------------

test_that("ee_dpsi_nu returns a square matrix of dimension K+1", {
  inp <- make_nu_inputs()
  out <- ee_dpsi_nu(inp$df, inp$nus)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(inp$K + 1L, inp$K + 1L))
})

test_that("ee_dpsi_nu is a diagonal matrix", {
  inp <- make_nu_inputs()
  out <- ee_dpsi_nu(inp$df, inp$nus)
  off_diag <- out[!diag(nrow(out))]
  expect_true(all(off_diag == 0))
})

test_that("ee_dpsi_nu diagonal entries equal -ns", {
  inp <- make_nu_inputs()
  out <- ee_dpsi_nu(inp$df, inp$nus)
  expect_equal(diag(out), rep(-inp$nus$ns, inp$K + 1L), tolerance = 1e-10)
})

test_that("ee_dpsi_nu is negative definite", {
  inp    <- make_nu_inputs()
  out    <- ee_dpsi_nu(inp$df, inp$nus)
  eigval <- eigen(out, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigval < 0))
})
