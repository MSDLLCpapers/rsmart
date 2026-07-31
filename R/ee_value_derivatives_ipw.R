#' Compute the derivative of the value estimating equation with respect to
#' propensity score parameters (IPW version)
#'
#' Computes \eqn{\partial \Psi_V / \partial \pi_k} for a single regime under
#' the inverse probability weighted (IPW) estimator. This is used in the An
#' matrix when no outcome regression models are specified.
#'
#' @param df A data frame containing the trial data, including treatment
#'   assignments, covariates, outcomes, and a \code{kappa} column.
#' @param pis A data frame of estimated propensity scores with columns
#'   \code{pi1, ..., piK}.
#' @param nus A list as returned by \code{\link{get_nu}}, containing the
#'   estimated stage probabilities.
#' @param regime_ind A matrix of regime consistency indicators for a single
#'   regime, with columns indicating whether each individual followed the
#'   regime at each stage.
#' @param p_fits A list of fitted propensity score model objects (one per stage).
#' @param k An integer indicating the stage index for which the derivative is
#'   computed.
#' @param qs Unused in the IPW version; included for interface consistency.
#'
#' @return A numeric row vector of partial derivatives of the value estimating
#'   equation with respect to the propensity score parameters at stage \code{k}.
#'
#' @noRd
ee_dpsiv_dpi_ipw <- function(df, pis, nus, regime_ind, p_fits, k, qs){
  K <- length(nus$nu) - 1
  
  Ck <- t(apply(regime_ind,1,cumprod))
  CR <- rowSums(Ck)  # consistency by coarsening of regime. 
  # if someone is coarsened due to CR, then their level should be odd
  R <- pmin(CR*2 + 1, df$kappa*2)  # get the minimum coarsening level for each individual
  # instead of infinity, we use levels K*2+1, (K+1)*2
  
  gprime <- function(xi, thetapi){
    # note that xi is 1x4 and the expfn is 1x1, result should be 4x1
    return (t(xi) %*% (exp(xi %*% thetapi) / ((1 + exp(xi %*% thetapi))^2)))
  }
  
  # 2K+1 level only
  # NOTE THAT BY COMPUTING FROM R=1 TO 2K, THE LAST DPSI.KR.DPIK IS ACTUALLY
  # FOR r=2K. WE USE THAT HERE
  r <- 2*K
  
  # now we get all the augmentation terms
  # the furthest stage reached with consistency is given as
  kr <- floor((r+1)/2)
  dfkm <- model.matrix(p_fits[[k]]@modelObj@model, data=df)
  # now consider the formation of the augmentation terms
  
  # note that lambdar is the pr(R=r|R \geq r)
  # K(r) = Pr(R > r) = Pr(R >= r+1)
  # We abuse the fact that I(R >=r) leads each augmentation term
  # So, Pr(R>=r) = K(r-1)
  # lambdar = [Pr(R=r) / Pr(R >= r)] = [Pr(R=r) / K(r-1)]
  # Pr(R=r) is given in paper supplemental
  # for odd r, this is the probability associated with NOT receiving regime 
  # for even r, this is the probability associated with kappa=k
  
  # note that individuals have pi=99 if they did not make it to a stage due
  # to time (kappa). Here, this can mean that some lambdar and Kr can take
  # strange values. However, the ifelse term should take care of that because
  
  # these are unchanged from the initial V estimation. 
  lambdar <-rep(1, nrow(df)) * (nus$nu[[kr]] - nus$nu[[kr+1]]) / nus$nu[[kr]] 
  Kr <- apply(pis[,1:kr, drop = FALSE], 1, prod) * nus$nu[[kr+1]]
  
  # here, if the stage k(r) is equal to the propensity of the stage considered, k,
  # we have dlambda is non-zero. 
  # else, it is 0. 
  dpsi.kr.dpik <- lapply(X=1L:nrow(df), FUN = function(i) (df[i,'kappa']>=kr) * 
                           (Kr[i] / pis[i,k]) * (2*Ck[i,k]-1) * (2*df[i,paste0('a',k)] -1) *
                           gprime(dfkm[i,,drop = FALSE], 
                                  p_fits[[k]]@fitObj$coefficients) ) 
  
  # NOTE THAT BY COMPUTING FROM R=1 TO 2K, THE LAST DPSI.KR.DPIK IS ACTUALLY
  # FOR r=2K. WE USE THAT HERE
  dInfTerm.dpik <- lapply(X=1L:nrow(df), FUN = function(i) - (dpsi.kr.dpik[[i]]) *
                            df[i,'y'] * (R[i] == (2*K + 1)) / (Kr[i]^2) ) 
  
  # we want to get dPsiV/dPi
  dV.dpik <- t(Reduce("+", dInfTerm.dpik))
  
  return (dV.dpik)
}

#' Compute the derivative of the value estimating equation with respect to
#' stage probability parameters (IPW version)
#'
#' Computes \eqn{\partial \Psi_V / \partial \nu_k} for a single regime under
#' the inverse probability weighted (IPW) estimator. This is used in the An
#' matrix when no outcome regression models are specified.
#'
#' @param df A data frame containing the trial data, including treatment
#'   assignments, covariates, outcomes, and a \code{kappa} column.
#' @param pis A data frame of estimated propensity scores with columns
#'   \code{pi1, ..., piK}.
#' @param nus A list as returned by \code{\link{get_nu}}, containing the
#'   estimated stage probabilities.
#' @param regime_ind A matrix of regime consistency indicators for a single
#'   regime.
#' @param k An integer indicating the nu index for which the derivative is
#'   computed.
#' @param qs Unused in the IPW version; included for interface consistency.
#'
#' @return A numeric scalar representing the partial derivative of the value
#'   estimating equation with respect to \eqn{\nu_k}.
#'
#' @noRd
ee_dpsiv_dnu_ipw <- function(df, pis, nus, regime_ind, k, qs){
  # note that nu has a single parameter for each k, 
  # since we compute this for each k, then this should be a scalar at the end
  
  K <- length(nus$nu) - 1
  
  Ck <- t(apply(regime_ind,1,cumprod))
  CR <- rowSums(Ck)  # consistency by coarsening of regime. 
  # if someone is coarsened due to CR, then their level should be odd
  R <- pmin(CR*2 + 1, df$kappa*2)  # get the minimum coarsening level for each individual
  # instead of infinity, we use levels K*2+1, (K+1)*2
  
  # iterate through the augmentation terms.
  # note nu1 only has presence in augmentation term 2, 
  # nu(K+1) appears in R=2K and 2K+1
  # all other nuk (k=2,...,K) appear in r=2k-2 (lambda,kr), 2k-1(Kr), 2k(lambda)
  
  # for terms r < 2k-1, the augmentation terms do not have pik, therefore, 
  # the derivative is 0. 
  
  dpsi.dnuk <- 0
  
  if (k == (K+1)){
    Kr <- apply(pis[,1:K, drop = FALSE], 1, prod) * nus$nu[[K+1]]
    dpsi.dnuk <- dpsi.dnuk + sum( - df[,'y'] * (R == (2*K + 1)) / (Kr * nus$nu[[K+1]] )  )
  }
  
  # we want to get dPsiV/dNu
  return (dpsi.dnuk)
}