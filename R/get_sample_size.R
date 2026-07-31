#' Determine sample size for a group sequential SMART design
#'
#' For specified operating characteristics, iteratively increases the sample
#' size until the desired power is achieved. Assumes the information fraction
#' remains unchanged as the sample size increases, which is reasonable for small
#' changes or at the design stage.
#'
#' @param variances A numeric vector of length \eqn{L} of variances of the 
#'   value estimators, i.e., \eqn{\sqrt{N} \times \mathrm{Cov}(\hat{\theta})}.
#'   These should reflect the population variances rather than sample variance
#'   or standard errors. 
#' @param beta A numeric value between 0 and 1 specifying the type II error
#'   rate. Power is \code{1 - beta}.
#' @param delta A numeric vector of length \eqn{L} of differences between
#'   regime values and the null value (or control arm).
#' @param bounds A numeric vector of length \eqn{S} with the stopping
#'   boundaries for analyses \eqn{1, \ldots, S}, or the boundaries from
#'   \code{\link{get_bounds}} function
#' @param n_init A positive integer specifying the initial total trial sample
#'   size \eqn{N} to begin the search.
#' @param corr A correlation matrix of dimension \eqn{L \times L} between
#'   regime value estimators across all analyses.
#' @param inf_frac A numeric vector of information fractions indicating when
#'   analyses are conducted.
#' @param n_split A numeric vector indicating the proportion of the sample
#'   size at the analysis times s=1,...,S. If no argument is given, assumes 
#'   the split is proportional to the information available. 
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
#'   }
#'
#' @export
get_sample_size <- function(variances, 
  beta, 
  delta, 
  bounds, 
  n_init = 100, 
  corr = bounds$corr, 
  inf_frac = bounds$inf_frac,
  n_split = NULL){
  
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

  ## Force evaluation of defaults before bounds is overwritten
  corr <- corr
  inf_frac <- inf_frac

  ## Validate corr and inf_frac after evaluation
  if (!is.matrix(corr) || nrow(corr) != ncol(corr)) {
    stop("corr must be a square matrix.")
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

  # stop if the correlation is not L by S
  if (S > 1){
    inf_frac.matrix <- sqrt(outer(inf_frac, inf_frac, pmin) / 
      outer(inf_frac, inf_frac, pmax))
    corr.inf <- kronecker(inf_frac.matrix, corr)
  } else {
    corr.inf <- corr
  }

  # this is used because of the possibility of interim analyses
  # having a sample size that is not proportional to information
  n <- ceiling(n_init * n_split)

  while (TRUE){
    mualt <- sqrt(c(rep(n, each=L) / variances)) * delta
    rejH0 <- 1 - pmvnorm(lower=-Inf, 
                         upper=rep(boundary, each = L), 
                         mean=mualt, 
                         corr = corr.inf)
    if (rejH0 < 1 - beta){
      n <- ceiling((max(n) + 1) * n_split)
    } else{
      break
    }
  }
  
  # calculate the expected cumulative rejections at each stage
  prop_rej <- c()
  for (ts in 1:S){
    exprej <- 1 - pmvnorm(lower=-Inf, 
                          upper=rep(boundary[1:ts], each=L),
                          mean = sqrt(c(rep(n[1:ts], each=L) / variances)) * delta,
                          corr=corr.inf[1:(L*ts), 1:(L*ts)])
    prop_rej <- c(prop_rej, exprej)
  }
  
  return(list("N" = n,
              "power" = rejH0,
              "prop_rej" = prop_rej,
              "variances" = variances,
              "beta" = beta,
              "delta" = delta,
              "bounds" = boundary,
              "n_init" = n_init,
              "corr" = corr,
              "inf_frac" = inf_frac,
              "n_split" = n_split))
}
