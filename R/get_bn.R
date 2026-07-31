#' Compute the Bn matrix for the sandwich variance estimator
#'
#' Computes the Bn matrix, defined as \eqn{n^{-1} \sum_{i=1}^{n} \Psi_i \Psi_i^T},
#' the empirical variance of the estimating equations. Each row of the
#' individual-level estimating equation matrix corresponds to one subject, with
#' columns ordered as \eqn{\pi_1, \ldots, \pi_K}, \eqn{\nu_1, \ldots, \nu_{K+1}},
#' \eqn{\beta^{(\ell)}_1, \ldots, \beta^{(\ell)}_K} for each regime
#' \eqn{\ell}, and the value estimating equations \eqn{V_1, \ldots, V_L}.
#' Based on Section 7 of Boos and Stefanski for the robust sandwich matrix.
#'
#' @param df A data frame containing the trial data, including treatment
#'   assignments, covariates, outcomes, and a \code{kappa} column.
#' @param p_fits A list of fitted propensity score model objects (one per stage),
#'   each a \code{modelObj} fit object.
#' @param nus A list as returned by \code{\link{get_nu}}, containing the
#'   estimated stage probabilities \code{nu}, the sample size \code{ns}, and
#'   \code{nd}.
#' @param q_all A list of fitted outcome regression objects for each regime,
#'   as returned by \code{\link{estimate_values}}. Can be an empty list when
#'   using IPW estimation.
#' @param dfs A list of data frames of value term components for each regime,
#'   as returned by \code{\link{estimate_values}}.
#' @param values A numeric vector of estimated regime values.
#' @param feasible_sets_indicator A logical value indicating whether feasible sets
#'   are present in the trial design (i.e., some treatments are deterministic
#'   based on response status).
#' @param q_list A list of outcome regression model specifications (one per
#'   stage), used to determine the number of Q-function models at each stage
#'   (e.g., separate models for responders and non-responders).
#'
#' @return A numeric matrix representing the Bn component of the sandwich
#'   variance estimator, with dimensions equal to the total number of estimated
#'   parameters.
#'
#' @export
get_bn <- function(df, p_fits, nus, q_all, dfs, values, feasible_sets_indicator, q_list){
  # each function returns psi_i as a row with the EE changing wrt parameters across the columns
  # therefore, the sum psi_i psi_i^T can be written instead as the matrix X^TX for X has rows 
  # xi = psi_i. We then divide the estimator by n-p following Boos Stefanski pg 320. 
  # if we divide by n, that is also asymptotically unbiased. 
  
  # get psi pi
  psi.pi <- ee_psi_pi(df, p_fits)
  
  # get psi nu
  psi.nu <- ee_psi_nu(df, nus)
  
  # get psi beta for all regimes 
  # the ee_psi_beta should do this for all of the embedded regimes 
  if (length(q_all) != 0){
    psi.beta <- ee_psi_beta(df, q_all, feasible_sets_indicator, q_list)
  } else{
    psi.beta <- NULL
  }
  
  # get values
  psi.v <- ee_psi_v(dfs, values, (df$kappa > 0))
  
  # psi <- as.data.frame(cbind(psi.pi1, psi.pi211, psi.pi201, psi.nu, psi.nu2, psi.nud, psibetas, psiv))
  psi <- as.data.frame(cbind(psi.pi, psi.nu, psi.beta, psi.v))
  
  # psi has dimension (n x p)
  
  # change made 3/7/2022: This should be n(s) not n(s) - p
  # return ( t(as.matrix(psi)) %*% as.matrix(psi) / (nus$ns - dim(psi)[2] ) )
  return ( t(as.matrix(psi)) %*% as.matrix(psi) / (nus$ns) )
  
}
