#' Compute the first chi-squared stopping boundary for a group sequential design
#'
#' Determines the first analysis stopping boundary based on a chi-squared
#' global test statistic. This is used when the comparison type is a global
#' chi-squared test (testing whether any regime differs from the others)
#' rather than individual regime comparisons against a fixed control.
#'
#' The function uses Monte Carlo simulation to estimate the boundary that
#' controls the familywise type I error rate at a specified level. A contrast
#' matrix \eqn{C} is constructed to form \eqn{L-1} linearly independent
#' comparisons among \eqn{L} regimes, and the chi-squared statistic is
#' computed as \eqn{(CZ)' (C \Sigma C')^{-1} (CZ)}.
#'
#' @param alpha A numeric value specifying the overall type I error rate to
#'   control. Default is 0.05.
#' @param inf_frac A numeric vector of information fractions indicating when
#'   analyses are conducted. Values should be between 0 and 1. The length
#'   determines the number of planned analyses \eqn{S}. Default is
#'   \code{c(0.5, 1)}.
#' @param spend_fn A character string specifying the alpha spending function.
#'   Currently used to adjust the \code{iota} scaling factors for each
#'   analysis boundary. For \code{"OF"} (O'Brien-Fleming), set \code{iota}
#'   to the information fractions. For \code{"Pocock"}, use
#'   \code{iota = rep(1, S)}. Default is \code{"OF"}.
#' @param corr A correlation matrix of the expected correlation between 
#'   the Z-statistics of regimes against a fixed value. Should have 
#'   dimension \eqn{L \times L} by \eqn{L \times L}, where \eqn{L} is
#'   the number of regimes. For \code{inf_frac} with length greater than one,
#'   the correlation structure across multiple analyses is computed. 
#' @param mu An optional numeric vector of means for the multivariate normal
#'   distribution for all regimes at a single analysis used in the Monte Carlo
#'   simulation. Default is \code{NULL}, which uses a zero vector (null
#'   hypothesis).
#' @param B A positive integer specifying the number of Monte Carlo samples.
#'   Default is 1000001.
#' @param seed An integer seed for reproducibility of the Monte Carlo
#'   simulation. Default is 1.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{bound}{A numeric vector of chi-squared stopping boundaries, one
#'       per analysis.}
#'     \item{spending}{A numeric vector of observed cumulative alpha spent at
#'       each analysis.}
#'     \item{typeI}{The overall achieved type I error rate across all
#'       analyses.}
#'     \item{convergence}{A logical value; always \code{TRUE} for Monte Carlo
#'       based computation.}
#'     \item{iters}{The number of Monte Carlo samples used (\code{B}).}
#'     \item{dfchi}{The degrees of freedom of the chi-squared statistic,
#'       equal to the rank of \eqn{C \Sigma C'}.}
#'     \item{alpha}{The input type I error rate.}
#'     \item{inf_frac}{The input information fractions.}
#'     \item{spend_fn}{The input spending function name.}
#'     \item{corr}{The input correlation matrix.}
#'     \item{mu}{The input mean vector.}
#'     \item{B}{The input number of Monte Carlo samples.}
#'     \item{seed}{The input random seed.}
#'   }
#'
#' @export
get_bounds_chi <- function(alpha = 0.05, inf_frac = c(0.5, 1),
                             spend_fn = "OF",
                             corr = diag(x = 1, nrow = 1, ncol = 1),
                             mu = NULL, B = 1000001, seed = 1) {

  # --- Input validation ---
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value in (0, 1).")
  }
  if (!is.numeric(inf_frac) || any(inf_frac <= 0) || any(inf_frac > 1) ||
      is.unsorted(inf_frac) || inf_frac[length(inf_frac)] != 1) {
    stop("inf_frac must be an increasing numeric vector in (0, 1] ending at 1.")
  }
  if (!is.character(spend_fn) || length(spend_fn) != 1 ||
      !spend_fn %in% c("OF", "Pocock")) {
    stop("spend_fn must be 'OF' or 'Pocock'.")
  }
  if (!is.matrix(corr) || nrow(corr) != ncol(corr)) {
    stop("corr must be a square matrix.")
  }
  if (nrow(corr) < 2) {
    stop("corr must be at least 2x2 (need >= 2 regimes for a chi-squared test).")
  }
  if (!is.null(mu)) {
    if (!is.numeric(mu) || length(mu) != nrow(corr)) {
      stop("mu must be a numeric vector with length equal to nrow(corr).")
    }
  }
  if (!is.numeric(B) || length(B) != 1 || B < 1 || B != round(B)) {
    stop("B must be a positive integer.")
  }
  if (!is.numeric(seed) || length(seed) != 1) {
    stop("seed must be a single numeric value.")
  }

  set.seed(seed)

  # Number of planned analyses
  S <- length(inf_frac)

  # Number of treatment regimes
  L <- dim(corr)[1] 

  # Create correlation matrix across analyses
  inf_frac.matrix <- sqrt(outer(inf_frac, inf_frac, pmin) / 
    outer(inf_frac, inf_frac, pmax))

  if (S > 1){
    corr.inf <- kronecker(inf_frac.matrix, corr)
  } else {
    corr.inf <- corr
  }

  # Determine the alpha spent at each of the analyses
  # OF: boundaries are scaled by information fraction (more conservative early)
  # Pocock: boundaries are constant across analyses (iota = 1)
  if (spend_fn == "OF"){
    alphat <- 2 - 2 * pnorm(qnorm(p=1-alpha/2, lower.tail=TRUE) / sqrt(inf_frac) )
  } else if (spend_fn == "Pocock") {
    alphat <- alpha * log(1 + (exp(1) -1 ) * inf_frac)
  }

  # Contrast matrix: L-1 independent comparisons (each regime vs last)
  Cmatrix <- cbind(diag(L - 1), -1)

  # Degrees of freedom = rank of C %*% Sigma %*% C'
  dfchi <- sum(svd(Cmatrix %*% corr.inf[1:L, 1:L] %*% t(Cmatrix))$d > 1e-8)

  # Default mean vector under the null hypothesis
  if (is.null(mu)) {
    mu <- rep(0, dim(corr.inf)[1])
  } else {
    mu <- rep(mu, S)
  }

  # --- Monte Carlo simulation ---
  # Generate B multivariate normal samples under mu
  Zb <- MASS::mvrnorm(n = B, mu = mu, Sigma = corr.inf)
  Zb <- matrix(Zb, nrow = B, byrow = FALSE)

  # Build block-diagonal contrast matrix for all S analyses
  CmatrixBdiag <- Cmatrix
  if (S > 1) {
    for (i in 2:S) {
      CmatrixBdiag <- Matrix::bdiag(CmatrixBdiag, Cmatrix)
    }
  }

  # Compute the block-diagonal covariance of contrasts: C %*% corr %*% C'
  bigCorr <- CmatrixBdiag %*% corr.inf %*% Matrix::t(CmatrixBdiag)
  bigCorrInv <- MASS::ginv(as.matrix(bigCorr))

  # Compute chi-squared statistics at each analysis for each MC sample
  # chiscores[b, s] = (C Z_s)' (C Sigma_s C')^{-1} (C Z_s)
  # column dums is faster than extracting the diagonal for large B
  chiscores <- matrix(0, nrow = B, ncol = S)
  for (i in 1:S) {
    row_idx <- ((i - 1) * L + 1):(i * L)
    contrast_idx <- ((i - 1) * (L - 1) + 1):(i * (L - 1))
    CZs <- Cmatrix %*% t(Zb)[row_idx, ]
    Sigma_s_inv <- MASS::ginv(as.matrix(bigCorr[contrast_idx, contrast_idx]))
    # chiscores[, i] <- diag(t(CZk) %*% Sigma_k_inv %*% CZk)
    chiscores[, i] <- colSums(CZs * (Sigma_s_inv %*% CZs))
  }

  # --- Cut scores at alpha spending ---
  # the boundaries at the cutoffs should result in the quantiles
  # the quantile of the first analysis should match the first alpha spend
  
  # now continue for the rest of the chiscores if there are more analyses
  # note that the global alpha spent should match the alphat at that stage
  # this means that we find the difference in the cumulative alpha spent 
  # and find the cutoff for the rejections at only that stage 
  alpha_inc <- diff(c(0, alphat))  # alpha expected rejections at each stage

  b1s <- unname(quantile(chiscores[, 1], probs = 1 - alpha_inc[1]))
  obs_spending <- mean(chiscores[, 1] > b1s)
  if (S > 1){
    # Track which samples have not yet been rejected
    not_rejected <- chiscores[, 1] < b1s[1]
    for (i in 2:S) {
      bi <- unname(quantile(chiscores[not_rejected, i], probs = 1 - alpha_inc[i]))
      b1s <- c(b1s, bi)
      # update observed alpha spending 
      obs_spending <- c(obs_spending, 
        obs_spending[i-1] + mean(chiscores[not_rejected, i] > bi) )
      # Update: also not rejected at this stage
      not_rejected <- not_rejected & (chiscores[, i] < bi)
    }
  }
  # Note that boundaries may not be exactly equal under Pocock-like bounds
  # this is not unique to Chi-square sequential testing

  # get observed type I for bounds
  confLevel <- 1 - mean(rowSums(chiscores > b1s) != 0) 

  return(list("bound" = b1s,
              "spending" = obs_spending, 
              "typeI" = 1 - confLevel,
              "convergence" = TRUE,
              "iters" = B,
              "dfchi" = dfchi,
              "alpha" = alpha,
              "inf_frac" = inf_frac,
              "spend_fn" = spend_fn,
              "corr" = corr,
              "mu" = mu,
              "B" = B,
              "seed" = seed))
}
