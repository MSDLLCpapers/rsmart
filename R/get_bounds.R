#' Compute all stopping boundaries for a two-analysis group sequential design
#'
#' Wrapper function that computes the stopping boundaries for both the first
#' and second analyses of a group sequential design with multiple treatment
#' regimes. Calls \code{\link{get_first_bound}} and \code{\link{get_next_bound}}
#' sequentially.
#'
#' @param alpha A numeric value specifying the overall type I error rate to
#'   control. Default is 0.05.
#' @param inf_frac A numeric vector of information fractions indicating when
#'   analyses are conducted.
#' @param spend_fn A character string specifying the alpha spending function.
#'   Either \code{"OF"} (O'Brien-Fleming) or \code{"Pocock"}.
#' @param corr A correlation matrix of Z-statistics at analysis time s. 
#'   Should have dimension \eqn{L \times L}.
#' @param test_type A character string specifying the type of test to be
#'   performed. Either \code{"one-sided"} or \code{"two-sided"}. Default is
#'   \code{"one-sided"}.
#' @param lambda A numeric value for the initial step size used in the
#'   iterative boundary search. Default is 0.1.
#' @param tol A numeric value specifying the convergence tolerance. Default is
#'   1e-6.
#' @param max_iter A positive integer specifying the maximum number of
#'   iterations. Default is 1000.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{bounds / bound}{A numeric value (single analysis) or numeric
#'       vector (sequential) of stopping boundaries. The key is \code{"bounds"}
#'       for a single analysis and \code{"bound"} for multiple analyses.}
#'     \item{spending}{A numeric value or vector of cumulative alpha spent at
#'       each analysis.}
#'     \item{convergence}{A logical value or vector indicating whether the
#'       algorithm converged at each analysis.}
#'     \item{iterations}{An integer value or vector of the number of
#'       iterations used at each analysis.}
#'     \item{alpha}{The input type I error rate.}
#'     \item{inf_frac}{The input information fractions.}
#'     \item{spend_fn}{The input spending function name.}
#'     \item{corr}{The input correlation matrix.}
#'     \item{test_type}{The input test type.}
#'   }
#'
#' @importFrom stats quantile
#' 
#' @export
get_bounds <- function(alpha = 0.05, inf_frac = c(0.5, 1), spend_fn = "OF",
                      corr = diag(x=1, nrow=1, ncol=1), test_type = "one-sided",
                      lambda = 0.1, tol = 1e-6, max_iter = 1000){
  
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
  if (!is.character(test_type) || length(test_type) != 1 ||
      !test_type %in% c("one-sided", "two-sided")) {
    stop("test_type must be 'one-sided' or 'two-sided'.")
  }
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0) {
    stop("lambda must be a positive numeric value.")
  }
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    stop("tol must be a positive numeric value.")
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter < 1 ||
      max_iter != round(max_iter)) {
    stop("max_iter must be a positive integer.")
  }

  # get the number of analyses and number of regimes
  S <- length(inf_frac)  # number of analyses to perform
  L <- dim(corr)[1]  # get number of regimes estimated each time point

  # iterate through the analyses
  b1s <- get_first_bound(alpha = alpha, inf_frac = inf_frac, spend_fn = spend_fn,
                       corr = corr, test_type = test_type,
                       lambda = lambda, tol = tol, max_iter = max_iter)
  
  # if information fraction is not length > 1, then we only need the first bound
  if (length(inf_frac) == 1) {
    return(list("bounds" = b1s$bound,
           "spending" = b1s$typeI,
           "convergence" = b1s$convergence,
           "iterations" = b1s$iters,
           "alpha" = alpha,
           "inf_frac" = inf_frac,
           "spend_fn" = spend_fn,
           "corr" = corr,
           "test_type" = test_type))
  }

  # Create correlation matrix across analyses
  inf_frac.matrix <- sqrt(outer(inf_frac, inf_frac, pmin) / 
    outer(inf_frac, inf_frac, pmax))

  if (S > 1){
    corr.inf <- kronecker(inf_frac.matrix, corr)
  } else {
    corr.inf <- corr
  }

  b2s <- get_next_bound(alpha = alpha, inf_frac = inf_frac, spend_fn = spend_fn,
                      corr = corr.inf, test_type = test_type,
                      lambda = lambda, tol = tol, 
                      prev_bound = rep(b1s$bound, L), 
                      s = 2, max_iter = max_iter)

  return(list("bound" = c(b1s$bound, b2s$bound), 
              "spending" = c(b1s$typeI, b2s$typeI),
              "convergence" = c(b1s$convergence, b2s$convergence),
              "iterations" = c(b1s$iters, b2s$iters),
              "alpha" = alpha,
              "inf_frac" = inf_frac,
              "spend_fn" = spend_fn,
              "corr" = corr,
              "test_type" = test_type) )
}
