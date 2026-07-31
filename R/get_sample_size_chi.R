
#' Determine sample size for a chi-squared global test in a group sequential
#' SMART design
#'
#' For specified operating characteristics, iteratively increases the sample
#' size until the desired power is achieved using a chi-squared global test
#' statistic. The chi-squared test assesses whether any treatment regime
#' differs from the others, rather than comparing individual regimes against
#' a fixed control.
#'
#' The function uses Monte Carlo simulation to estimate power at each
#' candidate sample size. A contrast matrix \eqn{C} is constructed to form
#' \eqn{L-1} linearly independent comparisons among \eqn{L} regimes, and
#' the chi-squared statistic is computed as
#' \eqn{(CZ)' (C \Sigma C')^{-1} (CZ)}.
#'
#' @param variances A numeric vector of length \eqn{L} of variances of the
#'   value estimators, i.e., \eqn{\sqrt{N} \times \mathrm{Cov}(\hat{\theta})}.
#'   These should reflect the population variances rather than sample variance
#'   or standard errors.
#' @param beta A numeric value between 0 and 1 specifying the type II error
#'   rate. Power is \code{1 - beta}.
#' @param delta A numeric vector of length \eqn{L} of change in the 
#'   regime values under the alternative hypothesis.
#' @param bounds A numeric vector of chi-squared stopping boundaries for
#'   analyses \eqn{1, \ldots, S}, or the output list from
#'   \code{\link{get_bounds_chi}}.
#' @param n_init A positive integer specifying the initial total trial sample
#'   size \eqn{N} to begin the search. Default is 100.
#' @param corr A correlation matrix of dimension \eqn{L \times L} between
#'   regime value estimators. Defaults to the correlation from the
#'   \code{bounds} list if provided.
#' @param inf_frac A numeric vector of information fractions indicating when
#'   analyses are conducted. Defaults to the information fractions from the
#'   \code{bounds} list if provided.
#' @param n_split A numeric vector indicating the proportion of the sample
#'   size at the analysis times \eqn{s = 1, \ldots, S}. If \code{NULL}
#'   (default), assumes the split is proportional to the information
#'   fractions.
#' @param seed An integer seed for reproducibility of the Monte Carlo
#'   simulation. Default is 2.
#' @param B A positive integer specifying the number of Monte Carlo samples
#'   for power estimation. Default is 10001.
#' @param lambda A positive numeric value specifying the incremental sample
#'   size increases when raising the sample size to find the required power 
#'   via simulation. After this is found, iteration from \eqn{n-lambda} to 
#'   \eqn{n+lambda} will be done to find the correct exact sample size. 
#' @param max_iter The maximum sample sizes n to evaluate. 
#'  
#' @return A list with the following components:
#'   \describe{
#'     \item{N}{A numeric vector of sample sizes at each analysis.}
#'     \item{power}{The achieved power at the final sample size.}
#'     \item{prop_rej}{A numeric vector of cumulative rejection probabilities
#'       at each analysis.}
#'     \item{variances}{The input variances.}
#'     \item{beta}{The input type II error rate.}
#'     \item{delta}{The input alternative differences.}
#'     \item{bounds}{The input stopping boundaries.}
#'     \item{n_init}{The input initial sample size.}
#'     \item{corr}{The input correlation matrix.}
#'     \item{inf_frac}{The input information fractions.}
#'     \item{n_split}{The input sample size split proportions.}
#'     \item{seed}{The input random seed.}
#'     \item{B}{The input number of Monte Carlo samples.}
#'     \item{lambda}{The input step size for the sample size search.}
#'   }
#'
#' @export
get_sample_size_chi <- function(
  variances, 
  beta, 
  delta, 
  bounds, 
  n_init = 100, 
  corr = bounds$corr, 
  inf_frac = bounds$inf_frac,
  n_split = NULL,
  seed = 2,
  B = 10001,
  lambda = 20,
  max_iter = 100){
  
  # --- Input validation ---
  if (!is.numeric(variances) || any(variances <= 0)) {
    stop("variances must be a positive numeric vector.")
  }
  if (!is.numeric(beta) || length(beta) != 1 || beta <= 0 || beta >= 1) {
    stop("beta must be a single numeric value in (0, 1).")
  }
  if (!is.numeric(delta) || any(!is.finite(delta))) {
    stop("delta must be a numeric vector with finite values.")
  }
  if (length(variances) != length(delta)) {
    stop("variances and delta must have the same length (number of regimes).")
  }
  if (!is.numeric(n_init) || length(n_init) != 1 || n_init < 1 ||
      n_init != round(n_init)) {
    stop("n_init must be a positive integer.")
  }
  if (!is.null(n_split)) {
    if (!is.numeric(n_split) || any(n_split <= 0) || any(n_split > 1)) {
      stop("n_split must be a numeric vector with values in (0, 1].")
    }
  }
  if (!is.numeric(seed) || length(seed) != 1) {
    stop("seed must be a single numeric value.")
  }
  if (!is.numeric(B) || length(B) != 1 || B < 1 || B != round(B)) {
    stop("B must be a positive integer.")
  }
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0) {
    stop("lambda must be a positive numeric value.")
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter < 1 ||
      max_iter != round(max_iter)) {
    stop("max_iter must be a positive integer.")
  }

  ## Force evaluation of defaults before bounds is overwritten
  corr <- corr
  inf_frac <- inf_frac

  ## Validate corr and inf_frac after evaluation
  if (!is.matrix(corr) || nrow(corr) != ncol(corr)) {
    stop("corr must be a square matrix.")
  }
  if (nrow(corr) < 2) {
    stop("corr must be at least 2x2 (need >= 2 regimes for a chi-squared test).")
  }
  if (!is.numeric(inf_frac) || any(inf_frac <= 0) || any(inf_frac > 1)) {
    stop("inf_frac must be a numeric vector with values in (0, 1].")
  }

  ## Extract boundary values from list input
  if (is.list(bounds)) {
    if (!is.null(bounds$bounds)) {
      boundary <- bounds$bounds
    } else if (!is.null(bounds$bound)) {
      boundary <- bounds$bound
    } else {
      stop("bounds list must contain a 'bounds' or 'bound' element.")
    }
  } else {
    boundary <- bounds
  }
  stopifnot(length(boundary) == length(inf_frac))


  if (is.null(n_split)) {
    n_split <- inf_frac
  }

  S <- length(boundary)

  L <- length(delta)

  # this is used because of the possibility of interim analyses
  # having a sample size that is not proportional to information
  n <- ceiling(n_init * n_split)
  
  # this is what changes as the sample changes
  # we use it to determine the alternative mean 
  mualt <- sqrt(c(rep(n, each=L) / variances)) * delta

  # Create correlation matrix across analyses
  inf_frac.matrix <- sqrt(outer(inf_frac, inf_frac, pmin) / 
    outer(inf_frac, inf_frac, pmax))

  if (S > 1){
    corr.inf <- kronecker(inf_frac.matrix, corr)
  } else {
    corr.inf <- corr
  }

  # create contrast matrices for the chi-squared tests and find 
  # associated degrees of freedom
  Cmatrix <- cbind(diag(L-1), -1) 
  dfchi <- sum(svd(Cmatrix %*% corr.inf[1:L,1:L] %*% Matrix::t(Cmatrix))$d > 1e-8)
  CmatrixBdiag <- Cmatrix
  if (S > 1){
    for (i in 2:S) {
      CmatrixBdiag <- bdiag(CmatrixBdiag, Cmatrix)
    }
  }

  getrej <- function(mu){
    ## now we want to get the chi squares with the new mu
    set.seed(seed)
    
    Zb <- MASS::mvrnorm(n=B, mu=mu, Sigma=corr.inf) # each row is a set of outcomes. 
    Zb <- matrix(Zb, nrow=B, byrow=FALSE)
    
    CZ <- CmatrixBdiag %*% t(Zb) # these are the contrasts
    bigCorr <- CmatrixBdiag %*% corr.inf %*% Matrix::t(CmatrixBdiag)
    bigCorrInv <- MASS::ginv(as.matrix(bigCorr))
    
    chiscores <- matrix(0, nrow=B, ncol=S)
    for (i in 1:S) {
      row_idx <- ((i - 1) * L + 1):(i * L)
      contrast_idx <- ((i - 1) * (L - 1) + 1):(i * (L - 1))
      CZs <- Cmatrix %*% t(Zb)[row_idx, ]
      Sigma_s_inv <- MASS::ginv(as.matrix(bigCorr[contrast_idx, contrast_idx]))
      # chiscores[, i] <- diag(t(CZk) %*% Sigma_k_inv %*% CZk)
      chiscores[, i] <- colSums(CZs * (Sigma_s_inv %*% CZs))
    }
    # the expected number of rejections is
    rejH0 <- mean(rowSums(chiscores > boundary) != 0)

    # get expected cumulative rejections at the stage
    prop_rej <- mean(chiscores[, 1] > boundary[1])
    if (S > 1){
      # Track which samples have not yet been rejected
      not_rejected <- chiscores[, 1] < boundary[1]
      for (i in 2:S) {
        prop_rej <- c(prop_rej, 
          prop_rej[i-1] + mean(chiscores[not_rejected, i] > boundary[i]) )
        # Update: also not rejected at this stage
        not_rejected <- not_rejected & (chiscores[, i] < boundary[i])
      }
    }

    return(list("rejH0" = rejH0, 
                "prop_rej" = prop_rej))
  }
  
  rejH0 <- getrej(mualt)
  
  if (rejH0$rejH0 > 1 - beta){
    print('At the initial sample size, the desired power is attained')
    return(list("N" = n,
            "power" = rejH0$rejH0,
            "prop_rej" = rejH0$prop_rej,
            "variances" = variances,
            "beta" = beta,
            "delta" = delta,
            "bounds" = boundary,
            "n_init" = n_init,
            "corr" = corr,
            "inf_frac" = inf_frac,
            "n_split" = n_split,
            "seed" = seed,
            "B" = B,
            "lambda" = lambda))
  }
  
  iters <- 1
  while(TRUE){
    while (TRUE){
      # in the case that the power is not enough yet, we want to increase the
      # sample size by lambda. Find new mean under alternative
      n <- ceiling((max(n) + lambda) * n_split)
      mualt <- sqrt(c(rep(n, each = L) / variances)) * delta
      # now we find power again
      rejH0 <- getrej(mualt)
      # check to see if we have achieved the power yet
      if (rejH0$rejH0 > 1 - beta){
        break  # exit the loop
      }
      iters <- iters + 1
      if (iters > max_iter) {
        message("Current power and sample size are: ", rejH0$rejH0, ", ", n)
        stop("Maximum iterations reached without achieving desired power.")
      }
    }
    # we reach this point when the power has been obtained.
    # now we want to go between the n2 -lambda and n2 by increments of 1
    n <- ceiling((max(n) - lambda) * n_split)
    while (TRUE){
      mualt <- sqrt(c(rep(n, each = L) / variances)) * delta
      # now we find power again
      rejH0 <- getrej(mualt)
      # once attained exit
      if (rejH0$rejH0 > 1 - beta){
        return(list("N" = n,
            "power" = rejH0$rejH0,
            "prop_rej" = rejH0$prop_rej,
            "variances" = variances,
            "beta" = beta,
            "delta" = delta,
            "bounds" = boundary,
            "n_init" = n_init,
            "corr" = corr,
            "inf_frac" = inf_frac,
            "n_split" = n_split,
            "seed" = seed,
            "B" = B,
            "lambda" = lambda))
      } else{
        # otherwise increase n2 by 1
        n <- ceiling((max(n) + 1) * n_split)
        iters <- iters + 1
        if (iters > max_iter) {
          message("Current power and sample size are: ", rejH0$rejH0, ", ", n)
          stop("Maximum iterations reached without achieving desired power.")
        }
      }
    }
  }
}