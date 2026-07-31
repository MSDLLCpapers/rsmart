#' Compute individual-level estimating equations for stage probabilities (nu)
#'
#' Computes \eqn{\Psi_{\nu,i}} for each individual \eqn{i}, the estimating
#' equation contributions for the stage arrival probabilities \eqn{\nu_k}. These
#' are used as part of the Bn matrix in the sandwich variance estimator.
#'
#' @param df A data frame containing the trial data, including a \code{kappa}
#'   column indicating the stage reached by each individual.
#' @param nus A list as returned by \code{\link{get_nu}}, containing the
#'   estimated stage probabilities \code{nu}, the sample size \code{ns}, and
#'   \code{nd}.
#'
#' @return A numeric matrix with rows corresponding to individuals and columns
#'   corresponding to each \eqn{\nu_k} parameter, representing the
#'   individual-level estimating equation contributions for nu.
#'
#' @noRd
ee_psi_nu <- function(df, nus){
  # psi nu for all k
  psi.nu <- c()
  for (k in (1:length(nus$nu))){
    psi.nu <- cbind(psi.nu, (df[,'kappa'] >= k) - (df[,'kappa'] >= 1)*nus$nu[[k]] )
  }
  # now for nud
  # psi.nu <- cbind(psi.nu, (df[,'kappa'] == (K+1)) - (df[,'kappa'] >= K)*nuFits$nd)
  return(psi.nu)  
}

#' Compute the negative derivative of the nu estimating equations with respect
#' to nu
#'
#' Computes \eqn{-\partial \Psi_{\nu} / \partial \nu}, the negative derivative
#' of the stage probability estimating equations with respect to the
#' \eqn{\nu_k} parameters. This forms a diagonal component of the An matrix in
#' the sandwich variance estimator.
#'
#' @param df A data frame containing the trial data, including a \code{kappa}
#'   column indicating the stage reached by each individual.
#' @param nus A list as returned by \code{\link{get_nu}}, containing the
#'   estimated stage probabilities \code{nu}, the sample size \code{ns}, and
#'   \code{nd}.
#'
#' @return A diagonal matrix with entries corresponding to
#'   \eqn{-\partial \Psi_{\nu_k} / \partial \nu_k} for each stage \eqn{k}.
#'
#' @noRd
ee_dpsi_nu <- function(df, nus){
  # begin with all k
  # this will now be block diagonal. 
  psi.nu <- c()
  for (k in (1:length(nus$nu))){
    psi.nu <- c(psi.nu, -sum(df[,'kappa'] >= 1) )
  }
  # now for nud
  # psi.nu <- cbind(psi.nu, (df[,'kappa'] >= K) )
  return(diag(psi.nu))  # the base function diag makes the vector a diagonal matrix. 
}
