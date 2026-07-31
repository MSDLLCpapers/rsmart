#' Compute stopping boundaries and sample size for a group sequential SMART
#'
#' Calculates the stopping boundaries and required sample size for a group
#' sequential design with multiple treatment regimes embedded in a SMART. The
#' boundaries are found first to control the type I error rate, then the sample
#' size is determined to achieve the desired power. This function wraps
#' \code{\link{get_bounds}} and \code{\link{get_sample_size}}.
#'
#' @param test_type A character string specifying the type of hypothesis test.
#'   A test type of either one-sided or two-sided. Two-sided testing must be
#'   symmetric. This input is ignored for a global Chi-squared test.
#' @param comp_type A character string specifying what comparison will be
#'   tested. Character string for either "fixed.control" or "global.chi.sq".
#' @param alpha A numeric value between 0 and 1 specifying the overall
#'   familywise type I error rate. Default is \code{0.05}.
#' @param beta A numeric value between 0 and 1 specifying the type II error
#'   rate. Power is \code{1 - beta}. Default is \code{0.1} for power of 0.9.
#' @param delta A numeric vector of length equal to the number of treatment
#'   regimes specifying the alternative differences between each regime's
#'   value and the null (or control) value.
#' @param n_init A positive integer specifying the initial total sample size
#'   to begin the iterative sample size search.
#' @param inf_frac A numeric vector of information fractions indicating when
#'   analyses are conducted. Values should be between 0 and 1, with the last
#'   element equal to 1. Default is \code{1} (single analysis). For example,
#'   \code{c(0.5, 1)} specifies an interim analysis at 50\\% information and a
#'   final analysis.
#' @param spend_fn A character string specifying the alpha spending function.
#'   Either \code{"OF"} (O'Brien-Fleming, default) or \code{"Pocock"}.
#' @param corr A correlation matrix of Z-statistics across all regimes and
#'   analyses. The dimension should be \eqn{L \times k} by \eqn{L \times k},
#'   where \eqn{L} is the number of treatment regimes. This matrix accounts
#'   for both within-analysis correlations between regimes and
#'   across-analysis correlations due to shared participants.
#' @param variances A numeric vector of length equal to the number of treatment
#'   regimes giving the variance of each value estimator scaled by sample
#'   size, i.e., \eqn{N \times \mathrm{Var}(\hat{V}_d)}.
#' @param lambdaB A numeric value for the initial step size used in the
#'   iterative boundary search. Default is 0.1.
#' @param lambdaC A numeric value for the initial step size used in the
#'   iterative sample size search. Default is 20.
#' @param tol A numeric value specifying the convergence tolerance. Default is
#'   1e-6.
#' @param max_iterB A positive integer specifying the maximum number of
#'   iterations for the boundary search. Default is 1000.
#' @param max_iterC A positive integer specifying the maximum number of
#'   iterations for the sample size search. Default is 100.
#' @param seed An integer seed for reproducibility of the Monte Carlo
#'   simulation. Default is 1.
#' @param Bb A positive integer specifying the number of Monte Carlo samples
#'   for the boundary computation. Default is 100001.
#' @param Bc A positive integer specifying the number of Monte Carlo samples
#'   for the sample size computation. Default is 10001.
#' @param mu A parameter vector for the multivariate normal distribution of
#'   the treatment effect mean, used only for the chi-square distribution
#'   in the case of a non-centrality assumption, though unlikely to be needed.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{boundaries}{The output from \code{\link{get_bounds}} (when
#'       \code{comp_type = "fixed.control"}) or \code{\link{get_bounds_chi}}
#'       (when \code{comp_type = "global.chi.sq"}), containing the stopping
#'       boundaries and associated metadata. See those functions for details.}
#'     \item{sample.size}{The output from \code{\link{get_sample_size}} or
#'       \code{\link{get_sample_size_chi}}, containing the required sample sizes,
#'       achieved power, cumulative rejection probabilities, and echoed input
#'       parameters. If \code{delta} is \code{NULL} (no alternative
#'       specified), the sample size is not computed and this element is
#'       \code{NULL}.}
#'     \item{inputs}{A list echoing all input parameters passed to
#'       \code{smart_design}.}
#'   }
#'
#' @export
smart_design <- function(
  test_type = "one-sided",
  comp_type = "fixed.control",
  alpha = 0.05,
  beta = 0.1, 
  delta = NULL,
  n_init = 400,
  inf_frac = 1,
  spend_fn = "OF",
  corr = NULL,
  variances = NULL,
  lambdaB = 0.1, 
  lambdaC = 20,
  tol = 1e-6,
  max_iterB = 1000,
  max_iterC = 100,
  seed = 1, 
  Bb = 100001,
  Bc = 10001,
  mu = NULL
){
  # --- Input validation ---
  if (!is.character(test_type) || length(test_type) != 1 ||
      !test_type %in% c("one-sided", "two-sided")) {
    stop("test_type must be 'one-sided' or 'two-sided'.")
  }
  if (!is.character(comp_type) || length(comp_type) != 1 ||
      !comp_type %in% c("fixed.control", "global.chi.sq")) {
    stop("comp_type must be 'fixed.control' or 'global.chi.sq'.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value in (0, 1).")
  }
  if (!is.numeric(beta) || length(beta) != 1 || beta <= 0 || beta >= 1) {
    stop("beta must be a single numeric value in (0, 1).")
  }
  if (!is.null(delta)) {
    if (!is.numeric(delta) || any(!is.finite(delta))) {
      stop("delta must be a numeric vector with finite values.")
    }
  }
  if (!is.numeric(n_init) || length(n_init) != 1 || n_init < 1 ||
      n_init != round(n_init)) {
    stop("n_init must be a positive integer.")
  }
  if (!is.numeric(inf_frac) || any(inf_frac <= 0) || any(inf_frac > 1) ||
      is.unsorted(inf_frac) || inf_frac[length(inf_frac)] != 1) {
    stop("inf_frac must be an increasing numeric vector in (0, 1] ending at 1.")
  }
  if (!is.character(spend_fn) || length(spend_fn) != 1 ||
      !spend_fn %in% c("OF", "Pocock")) {
    stop("spend_fn must be 'OF' or 'Pocock'.")
  }
  if (is.null(corr)) {
    stop("corr must be provided as a square correlation matrix.")
  }
  if (!is.matrix(corr) || nrow(corr) != ncol(corr)) {
    stop("corr must be a square matrix.")
  }
  if (!is.null(variances)) {
    if (!is.numeric(variances) || any(variances <= 0)) {
      stop("variances must be a positive numeric vector.")
    }
  }
  if (!is.null(delta) && is.null(variances)) {
    stop("variances must be provided when delta is specified.")
  }
  if (!is.numeric(lambdaB) || length(lambdaB) != 1 || lambdaB <= 0) {
    stop("lambdaB must be a positive numeric value.")
  }
  if (!is.numeric(lambdaC) || length(lambdaC) != 1 || lambdaC <= 0) {
    stop("lambdaC must be a positive numeric value.")
  }
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    stop("tol must be a positive numeric value.")
  }
  if (!is.numeric(max_iterB) || length(max_iterB) != 1 || max_iterB < 1 ||
      max_iterB != round(max_iterB)) {
    stop("max_iterB must be a positive integer.")
  }
  if (!is.numeric(max_iterC) || length(max_iterC) != 1 || max_iterC < 1 ||
      max_iterC != round(max_iterC)) {
    stop("max_iterC must be a positive integer.")
  }
  if (!is.numeric(seed) || length(seed) != 1) {
    stop("seed must be a single numeric value.")
  }
  if (!is.numeric(Bb) || length(Bb) != 1 || Bb < 1 || Bb != round(Bb)) {
    stop("Bb must be a positive integer.")
  }
  if (!is.numeric(Bc) || length(Bc) != 1 || Bc < 1 || Bc != round(Bc)) {
    stop("Bc must be a positive integer.")
  }

  # note that for any regime against a fixed control, the boundaries can be 
  # calculated using get_bounds() to sequentially estimate them. 

  # boundaries are iteratively calculated using the starting point of a single
  # analysis for a single treatment regime. Then expanded for multiple analyses
  # using the alpha spent at that time and the correlation structure between 
  # groups. 

  if (comp_type == "fixed.control"){
    # for fixed control, we can use the boundaries directly
    boundaries <- get_bounds(alpha = alpha, inf_frac = inf_frac,
      spend_fn = spend_fn, corr = corr, test_type = test_type,
      lambda = lambdaB, tol = tol, max_iter = max_iterB)
    if (!is.null(delta)){
      sample.size <- get_sample_size(variances = variances, 
                beta = beta, 
                delta = delta, 
                bounds = boundaries, 
                n_init = n_init, 
                n_split = NULL)
    } else{
      sample.size = NULL
    }
  } else if (comp_type == "global.chi.sq") {
    # for global chi-squared, we need to adjust the boundaries
    boundaries <- get_bounds_chi(alpha = alpha, inf_frac = inf_frac,
      spend_fn = spend_fn, corr = corr, mu = mu, B = Bb, seed = seed)
    if (!is.null(delta)){
      sample.size <- get_sample_size_chi(variances = variances, 
                beta = beta, 
                delta = delta, 
                bounds = boundaries, 
                n_init = n_init, 
                n_split = NULL,
                seed = seed,
                B = Bc,
                lambda = lambdaC,
                max_iter = max_iterC)
    } else{
      sample.size = NULL
    }
  }

  # create the return information
  return(list("boundaries" = boundaries, 
  "sample.size" = sample.size,
  "inputs" = list(
    "test_type" = test_type,
    "comp_type" = comp_type,
    "alpha" = alpha,
    "beta" = beta,
    "delta" = delta,
    "n_init" = n_init,
    "inf_frac" = inf_frac,
    "spend_fn" = spend_fn,
    "corr" = corr,
    "variances" = variances,
    "lambdaB" = lambdaB,
    "lambdaC" = lambdaC,
    "tol" = tol,
    "max_iterB" = max_iterB,
    "max_iterC" = max_iterC,
    "seed" = seed,
    "Bb" = Bb,
    "Bc" = Bc,
    "mu" = mu
  )))
}

