#' Assemble the full derivative matrix of value estimating equations (IPW
#' version)
#'
#' Constructs the matrix \eqn{-\partial \Psi_V / \partial \theta} for the
#' inverse probability weighted (IPW) estimator across all regimes. Each row
#' corresponds to one regime's value estimating equation derivative, with
#' columns ordered as propensity score parameters, nu parameters, and value
#' parameters. This is used in the An matrix when no outcome regression models
#' are specified.
#'
#' @param df A data frame containing the trial data, including treatment
#'   assignments, covariates, outcomes, and a \code{kappa} column.
#' @param pis A data frame of estimated propensity scores with columns
#'   \code{pi1, ..., piK}.
#' @param nus A list as returned by \code{\link{get_nu}}, containing the
#'   estimated stage probabilities.
#' @param regime_all A list of regime objects, each containing a \code{regime}
#'   matrix and a \code{regime_ind} indicator matrix.
#' @param p_fits A list of fitted propensity score model objects (one per stage).
#' @param q_all A list of fitted outcome regression objects (expected to be
#'   empty for IPW); included for interface consistency.
#' @param feasible_sets_indicator A logical value indicating whether feasible sets
#'   are present in the trial design.
#'
#' @return A numeric matrix with rows corresponding to regimes and columns
#'   corresponding to all estimated parameters, representing
#'   \eqn{-\partial \Psi_V / \partial \theta} under the IPW estimator.
#'
#' @noRd
ee_dpsiv_ipw <- function(df, pis, nus, regime_all, p_fits, q_all, feasible_sets_indicator){
  # note that we need to get a row vector for each ell with the 
  # row having entries in the same order as Bn
  # i.e. Pi1-PiK, nu1-nuK+1, ell1 Beta1-BetaK, ellL Beta1-BetaK, V1-Vell. 
  # for ell != current, the betas should be zero. 
  # since V1-Vell is done as a block matrix, we can include this at the very end. 
  
  for (ell in 1:length(regime_all)){
    # start with pi
    dPsiV.dPi.ell <- c()
    dPsiV.dNu.ell <- c()
    for (k in 1:length(p_fits)){
      dPsiV.dPi.ell <- cbind(dPsiV.dPi.ell, 
                             matrix(
                               ee_dpsiv_dpi_ipw(df=df, pis=pis, nus=nus, 
                                               regime_ind=regime_all[[ell]]$regime_ind, 
                                               p_fits=p_fits, k=k, qs=NULL) 
                               ,nrow=1)
      )
      
      dPsiV.dNu.ell <- cbind(dPsiV.dNu.ell, 
                             matrix(
                               ee_dpsiv_dnu_ipw(df=df, pis=pis, nus=nus, 
                                               regime_ind=regime_all[[ell]]$regime_ind,
                                               k=k, qs=NULL) 
                               ,nrow=1)
      )
      if (k==length(p_fits)){
        dPsiV.dNu.ell <- cbind(dPsiV.dNu.ell, 
                               matrix(
                                 ee_dpsiv_dnu_ipw(df=df, pis=pis, nus=nus, 
                                                 regime_ind=regime_all[[ell]]$regime_ind,
                                                 k=k+1, qs=NULL) 
                                 ,nrow=1)
        )
      }
      
    }# end for loop over k stages
    if (ell==1){
      # initialize dPsiV as just that of Pi, Nu, Beta for regime 1
      dPsiV <- cbind(dPsiV.dPi.ell, dPsiV.dNu.ell)  
    } else{
      dPsiV.ell <- cbind(dPsiV.dPi.ell, dPsiV.dNu.ell)
      # now the columns line up to bind the rows together
      dPsiV <- rbind(dPsiV, dPsiV.ell)
    }
    
  }# end for loop over ell regime number 
  # now we should have the matrix of dPsiV / dPi dNu dBeta 11 - dBetaLK. 
  # add in the diagonal matrix corresponding to V
  dPsiV.V <- ee_dpsiv_dv(nregimes=length(regime_all), nus=nus)
  
  dPsiV <- cbind(dPsiV, dPsiV.V)
  
  # this should return an entire row. It is still -dPsiV/dTheta
  # It is a row rather than a block because the estimating equations of V
  # are functions of the other parameters
  return (dPsiV) 
  
  
}
