#' Compute the An matrix for the sandwich variance estimator
#'
#' Computes the An matrix, defined as \eqn{-\partial \Psi / \partial \theta},
#' where \eqn{\Psi} is the stacked estimating equations vector. The estimating
#' equations are ordered to match Bn: \eqn{\pi_1, \ldots, \pi_K},
#' \eqn{\nu_1, \ldots, \nu_{K+1}}, \eqn{\beta^{(\ell)}_1, \ldots, \beta^{(\ell)}_K}
#' for each regime \eqn{\ell}, and the value estimating equations
#' \eqn{V_1, \ldots, V_L}. Based on Section 7 of Boos and Stefanski for the
#' robust sandwich matrix.
#'
#' @param df A data frame containing the trial data, including treatment
#'   assignments, covariates, outcomes, and a \code{kappa} column.
#' @param pis A data frame of estimated propensity scores with columns
#'   \code{pi1, ..., piK}, one column per stage.
#' @param p_fits A list of fitted propensity score model objects (one per stage),
#'   each a \code{modelObj} fit object.
#' @param nus A list as returned by \code{\link{get_nu}}, containing the
#'   estimated stage probabilities \code{nu}, the sample size \code{ns}, and
#'   \code{nd}.
#' @param q_all A list of fitted outcome regression objects for each regime,
#'   as returned by \code{\link{estimate_values}}. Can be an empty list when
#'   using IPW estimation.
#' @param values A numeric vector of estimated regime values.
#' @param regime_all A list of regime objects, each containing a \code{regime}
#'   matrix and a \code{regime_ind} indicator matrix.
#' @param dfs A list of data frames of value term components for each regime,
#'   as returned by \code{\link{estimate_values}}.
#' @param feasible_sets_indicator A logical value indicating whether feasible sets
#'   are present in the trial design (i.e., some treatments are deterministic
#'   based on response status).
#' @param q_list A list of outcome regression model specifications (one per
#'   stage), used to determine the number of Q-function models at each stage
#'   (e.g., separate models for responders and non-responders).
#'
#' @return A numeric matrix representing the An component of the sandwich
#'   variance estimator, with dimensions equal to the total number of estimated
#'   parameters.
#'
#' @export
get_an <- function(df, pis, p_fits, nus, q_all, values, regime_all, dfs, 
  feasible_sets_indicator, q_list){
  
  # for pi, nu, and beta, these are the -derivatives of the estimating equations 
  # of those parameters with respect to only those parameters. This is because
  # they are 0 for all other entries.
  # however, dpsi.v is a matrix with each row as the -derivative of the regime-
  # specific estimating equation with respect to all parameters. 
  # so when we combine these terms, we create a block diagonal matrix 
  # from the first three (pi, nu, beta), then bind the rows of v to the bottom. 
  # we must bind a matrix of zeros to the columns of the (pi, nu, beta) because
  # it does not have the columns for dpsi.v
  
  # get dpsi pi, one call for all k and ell; this is NOT -dpsi.pi/dpi
  dpsi.pi <- ee_dpsi_pi(df, p_fits)
  
  # get dpsi nu, one call for all k and ell; this is NOT -dpsi.nu/dnu
  dpsi.nu <- ee_dpsi_nu(df, nus)
  
  # get psi beta for all regimes 
  # the ee_psi_beta should do this for all of the embedded regimes 
  # the regime_all is added because now we deal with different design matrices
  # this is NOT -dpsi.beta/dbeta
  if (length(q_all) != 0){
    dpsi.beta <- ee_dpsi_beta(df, q_all, regime_all, feasible_sets_indicator, q_list)
    # get values ee derivatives
    dpsi.v <- ee_dpsiv(df, pis, nus, regime_all, p_fits, q_all, feasible_sets_indicator, q_list)
  } else{
    dpsi.v <- ee_dpsiv_ipw(df, pis, nus, regime_all, p_fits, q_all, feasible_sets_indicator)
  }
  
  # create block diagonal matrix 
  if (length(q_all) != 0){
    dpsi.pi.nu.beta <- Matrix::bdiag(dpsi.pi, dpsi.nu, dpsi.beta)
  } else{
    dpsi.pi.nu.beta <- Matrix::bdiag(dpsi.pi, dpsi.nu)
  }
  
  # divide matrix by the sample size nus$ns, 
  return (rbind(cbind(dpsi.pi.nu.beta,
              matrix(0, nrow=nrow(dpsi.pi.nu.beta), ncol=nrow(dpsi.v)) 
              ),
              dpsi.v) / nus$ns
  )
  
}
