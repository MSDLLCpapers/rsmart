#' Compute individual-level estimating equations for value parameters
#'
#' Computes \eqn{\Psi_{V,i}} for each individual \eqn{i}, the estimating
#' equation contributions for the regime value parameters. These are used as
#' part of the Bn matrix in the sandwich variance estimator.
#'
#' @param dfs A list of matrices of value term components for each regime, as
#'   returned by \code{\link{estimate_values}}.
#' @param values A numeric vector of estimated regime values.
#' @param cKappa A logical vector indicating whether each individual has been
#'   observed (i.e., \code{kappa > 0}).
#'
#' @return A numeric matrix with rows corresponding to individuals and columns
#'   corresponding to each regime value parameter, representing the
#'   individual-level estimating equation contributions for V.
#'
#' @noRd
ee_psi_v <- function(dfs, values, cKappa){
  # cKappa indicates if individual observed yet (df$kappa > 0)
  makeNice <- function(...) {matrix(unlist(...), nrow=length(...), byrow=TRUE)}
  return ( t( makeNice(lapply(X=1:length(dfs), FUN=function(i) rowSums(dfs[[i]]) - values[i]*cKappa)
  )
  )
  )
}

#' Assemble the full derivative matrix of value estimating equations (AIPW
#' version)
#'
#' Constructs the matrix \eqn{-\partial \Psi_V / \partial \theta} for the
#' augmented inverse probability weighted (AIPW) estimator across all regimes.
#' Each row corresponds to one regime's value estimating equation derivative,
#' with columns ordered as propensity score parameters, nu parameters, beta
#' parameters (by regime), and value parameters.
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
#' @param q_all A list of fitted outcome regression objects for each regime,
#'   as returned by \code{\link{estimate_values}}.
#' @param feasible_sets_indicator A logical value indicating whether feasible sets
#'   are present in the trial design.
#' @param q_list A list of outcome regression model specifications (one per
#'   stage), used to determine the number of Q-function models at each stage
#'   (e.g., separate models for responders and non-responders).
#'
#' @return A numeric matrix with rows corresponding to regimes and columns
#'   corresponding to all estimated parameters, representing
#'   \eqn{-\partial \Psi_V / \partial \theta} under the AIPW estimator.
#'
#' @noRd
ee_dpsiv <- function(df, pis, nus, regime_all, p_fits, q_all, feasible_sets_indicator, q_list){
  # note that we need to get a row vector for each ell with the 
  # row having entries in the same order as Bn
  # i.e. Pi1-PiK, nu1-nuK+1, ell1 Beta1-BetaK, ellL Beta1-BetaK, V1-Vell. 
  # for ell != current, the betas should be zero. 
  # since V1-Vell is done as a block matrix, we can include this at the very end. 
  
  for (ell in 1:length(regime_all)){
    # start with pi
    dPsiV.dPi.ell <- c()
    dPsiV.dNu.ell <- c()
    dPsiV.dBeta.ell <- c()
    for (k in 1:length(p_fits)){
      dPsiV.dPi.ell <- cbind(dPsiV.dPi.ell, 
                             matrix(
                             ee_dpsiv_dpi(df=df, pis=pis, nus=nus, 
                                         regime_ind=regime_all[[ell]]$regime_ind, 
                                         p_fits=p_fits, k=k, qs=q_all[[ell]]) 
                             ,nrow=1)
      )
    
      dPsiV.dNu.ell <- cbind(dPsiV.dNu.ell, 
                             matrix(
                             ee_dpsiv_dnu(df=df, pis=pis, nus=nus, 
                                         regime_ind=regime_all[[ell]]$regime_ind,
                                         k=k, qs=q_all[[ell]]) 
                             ,nrow=1)
      )
      if (k==length(p_fits)){
        dPsiV.dNu.ell <- cbind(dPsiV.dNu.ell, 
                               matrix(
                                 ee_dpsiv_dnu(df=df, pis=pis, nus=nus, 
                                             regime_ind=regime_all[[ell]]$regime_ind,
                                             k=k+1, qs=q_all[[ell]]) 
                                 ,nrow=1)
        )
      }
      dPsiV.dBeta.ell <- cbind(dPsiV.dBeta.ell,
                               matrix(
                               ee_dpsiv_dbeta(df=df, pis=pis, nus=nus, 
                                             regime_ind=regime_all[[ell]]$regime_ind,
                                             regime = regime_all[[ell]]$regime,
                                             q_fits = q_all[[ell]]$q_fits,
                                             k=k, feasible_sets_indicator=feasible_sets_indicator,
                                             q_list = q_list) 
                               ,nrow=1)
      )
                    
    }# end for loop over k stages
    if (ell==1){
      # initialize dPsiV as just that of Pi, Nu, Beta for regime 1
      dPsiV <- cbind(dPsiV.dPi.ell, dPsiV.dNu.ell, dPsiV.dBeta.ell)  
      zeroforbeta <- ncol(dPsiV.dBeta.ell)
    } else{
      # add a columns of zeros equal to the dimension of new betas to dPsiV
      dPsiV <- cbind(dPsiV, 
                     matrix(0, nrow=nrow(dPsiV), ncol=ncol(dPsiV.dBeta.ell)))
      dPsiV.ell <- cbind(dPsiV.dPi.ell, dPsiV.dNu.ell, matrix(0, ncol=zeroforbeta), 
                         dPsiV.dBeta.ell)
      # now the columns line up to bind the rows together
      dPsiV <- rbind(dPsiV, dPsiV.ell)
      # update the count of betas from previous regimes 
      zeroforbeta <- zeroforbeta + ncol(dPsiV.dBeta.ell)
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
