# This script defines the functions to return the partial derivatives of the 
# terms associated with - dPsiV / dTheta. Each of these should return the 
# vector of associated derivatives for that function

#### With respect to pi ####
# note that pik only exists in terms for which kr >= k, so we need only iterate
# from the first r such that kr >=k, which is r=2k-1 and on

#' Compute the derivative of the value estimating equation with respect to
#' propensity score parameters (AIPW version)
#'
#' Computes \eqn{\partial \Psi_V / \partial \pi_k} for a single regime under
#' the augmented inverse probability weighted (AIPW) estimator. Iterates
#' through the augmentation terms from \eqn{r = 2k-1} to \eqn{2K} and the
#' IPW (infinity) term.
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
#' @param qs A list containing outcome regression fitted values for a single
#'   regime, as returned by \code{\link{get_q_fits}}.
#'
#' @return A numeric row vector of partial derivatives of the value estimating
#'   equation with respect to the propensity score parameters at stage \code{k}.
#'
#' @noRd
ee_dpsiv_dpi <- function(df, pis, nus, regime_ind, p_fits, k, qs){
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
  
  # iterate through the augmentation terms.
  # for terms r < 2k-1, the augmentation terms do not have pik, therefore, 
  # the derivative is 0. 
  dAugTerms.dpik <- lapply(vector(mode="list", length=nrow(df)), function(x) 0)
  for (r in (2*k-1):(2*K)){
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
    if (r %% 2 == 1){
      # odd r
      lambdar <- ( (pis[,kr]*(1-Ck[,kr]) + (1-pis[,kr])*Ck[,kr]) )  
      if (r==1){
        Kr <- ( (pis[,kr]*Ck[,kr] + (1-pis[,kr])*(1-Ck[,kr])) )  
      } else {
        Kr <- apply(pis[,1:(kr-1), drop = FALSE], 1, prod) * nus$nu[[kr]] * 
          ( (pis[,kr]*Ck[,kr] + (1-pis[,kr])*(1-Ck[,kr])) )  
      }
    } else{
      # even r
      # this is 0 for AIPW estimator 
      # lambdar needs to be a vector with the same number of rows of df
      lambdar <-rep(1, nrow(df)) * (nus$nu[[kr]] - nus$nu[[kr+1]]) / nus$nu[[kr]] 
      Kr <- apply(pis[,1:kr, drop = FALSE], 1, prod) * nus$nu[[kr+1]]
    }
    
    # here, if the stage k(r) is equal to the propensity of the stage considered, k,
    # we have dlambda is non-zero. 
    # else, it is 0. 
    if ((kr == k) & (r %% 2 == 1)){
      dpsi.lambdar.dpik <- lapply(X=1L:nrow(df), FUN = function(i) (df[i,'kappa']>=kr) * 
                                    (1-2*Ck[i,kr]) * (2*df[i,paste0('a',kr)] -1) *
                                    gprime(dfkm[i,,drop = FALSE], 
                                           p_fits[[k]]@fitObj$coefficients) 
      )
    } else{
      # create a list of all 0s so that the lapply later will still work. 
      dpsi.lambdar.dpik <- lapply(vector(mode="list", length=nrow(df)), function(x) 0)
    }
    if (r == 1){
      dpsi.kr.dpik <- lapply(X=1L:nrow(df), FUN = function(i) (df[i,'kappa']>=kr) * 
                               (2*Ck[i,k]-1) * (2*df[i,paste0('a',k)] -1) *
                               gprime(dfkm[i,,drop = FALSE], 
                                      p_fits[[k]]@fitObj$coefficients) ) 
    } else{
      dpsi.kr.dpik <- lapply(X=1L:nrow(df), FUN = function(i) (df[i,'kappa']>=kr) * 
                               (Kr[i] / pis[i,k]) * (2*Ck[i,k]-1) * (2*df[i,paste0('a',k)] -1) *
                               gprime(dfkm[i,,drop = FALSE], 
                                      p_fits[[k]]@fitObj$coefficients) ) 
    } 
    # this is only the rth augmentation term
    dOneAugTerms.dpik <- lapply(X=1L:nrow(df), FUN = function(i)  ( (Kr[i] * (-dpsi.lambdar.dpik[[i]] * (R[i] >= r) ) - 
                                                                       ((R[i] == r) - (R[i] >= r)*lambdar[i]) * dpsi.kr.dpik[[i]] ) / 
                                                                      (Kr[i]^2) ) *
                                  qs$mod_regime_vhats[i,kr] )
    # add the rth augmentation term to the list
    dAugTerms.dpik <- lapply(X=1L:nrow(df), FUN = function(i) dAugTerms.dpik[[i]] + dOneAugTerms.dpik[[i]])
  } # end for loop  
  # NOTE THAT BY COMPUTING FROM R=1 TO 2K, THE LAST DPSI.KR.DPIK IS ACTUALLY
  # FOR r=2K. WE USE THAT HERE
  dInfTerm.dpik <- lapply(X=1L:nrow(df), FUN = function(i) - (dpsi.kr.dpik[[i]]) *
                            df[i,'y'] * (R[i] == (2*K + 1)) / (Kr[i]^2) ) 
  
  # we want to get dPsiV/dPi
  dV.dpik <- t(Reduce("+", dAugTerms.dpik) + Reduce("+", dInfTerm.dpik))
  
  return (dV.dpik)
}


#### With respect to nu ####
# k is meant to index nu.

#' Compute the derivative of the value estimating equation with respect to
#' stage probability parameters (AIPW version)
#'
#' Computes \eqn{\partial \Psi_V / \partial \nu_k} for a single regime under
#' the augmented inverse probability weighted (AIPW) estimator. Iterates
#' through all augmentation terms and the IPW (infinity) term.
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
#' @param qs A list containing outcome regression fitted values for a single
#'   regime, as returned by \code{\link{get_q_fits}}.
#'
#' @return A numeric scalar representing the partial derivative of the value
#'   estimating equation with respect to \eqn{\nu_k}.
#'
#' @noRd
ee_dpsiv_dnu <- function(df, pis, nus, regime_ind, k, qs){
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
  
  for (r in 1:(2*K)){
    # now we get all the augmentation terms
    # the furthest stage reached with consistency is given as
    kr <- floor((r+1)/2)
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
    if (r %% 2 == 1){
      # odd r
      lambdar <- ( (pis[,kr]*(1-Ck[,kr]) + (1-pis[,kr])*Ck[,kr]) )  
      if (r==1){
        Kr <- ( (pis[,kr]*Ck[,kr] + (1-pis[,kr])*(1-Ck[,kr])) )  
      } else {
        Kr <- apply(pis[,1:(kr-1), drop = FALSE], 1, prod) * nus$nu[[kr]] * 
          ( (pis[,kr]*Ck[,kr] + (1-pis[,kr])*(1-Ck[,kr])) )  
      }
      if ((kr == k) & (k >1)){
        dpsi.kr.nuk <- Kr / nus$nu[[k]]  
      } else{
        dpsi.kr.nuk <- 0
      }
      dpsi.lambdar.nuk <- 0
    } else{
      # even r
      # this is 0 for AIPW estimator 
      lambdar <-rep(1, nrow(df)) * (nus$nu[[kr]] - nus$nu[[kr+1]]) / nus$nu[[kr]] 
      Kr <- apply(pis[,1:kr, drop = FALSE], 1, prod) * nus$nu[[kr+1]]
      if(kr==k){
        dpsi.lambdar.nuk <- nus$nu[[k+1]] * (nus$nu[[k]]^(-2))
        dpsi.kr.nuk <- 0
      }  
      if (kr==(k-1)){
        dpsi.lambdar.nuk <- -(nus$nu[[kr]]^(-1)) 
        dpsi.kr.nuk <- Kr / nus$nu[[kr+1]] 
      }
      
    }
    # do the augmentation term here...
    dpsi.dnuk <- dpsi.dnuk + sum(
      qs$mod_regime_vhats[,kr] * 
        ( (Kr * (-dpsi.lambdar.nuk) * (R >= r)) - 
            (( (R == r) - (lambdar*(R >= r)) ) * dpsi.kr.nuk) ) / (Kr^2), na.rm=TRUE)
    
  } # end for loop  
  
  if (k == (K+1)){
    Kr <- apply(pis[,1:K, drop = FALSE], 1, prod) * nus$nu[[K+1]]
    dpsi.dnuk <- dpsi.dnuk + sum( - df[,'y'] * (R == (2*K + 1)) / (Kr * nus$nu[[K+1]] )  )
  }
  
  # we want to get dPsiV/dNu
  return (dpsi.dnuk)
}

#### With respect to beta ####

# for the Q functions, these are only going to be non-zero for when k(r) = k,
# where k is the index of the q function considered.
# These differ by ell because it will use the modified design matrix.

#' Compute the derivative of the value estimating equation with respect to
#' outcome regression coefficients (beta)
#'
#' Computes \eqn{\partial \Psi_V / \partial \beta_k} for a single regime under
#' the augmented inverse probability weighted (AIPW) estimator. The Q-function
#' coefficients at stage \eqn{k} only appear in the augmentation terms
#' \eqn{r = 2k-1} and \eqn{r = 2k}.
#'
#' @param df A data frame containing the trial data, including treatment
#'   assignments, covariates, outcomes, and a \code{kappa} column.
#' @param pis A data frame of estimated propensity scores with columns
#'   \code{pi1, ..., piK}.
#' @param nus A list as returned by \code{\link{get_nu}}, containing the
#'   estimated stage probabilities.
#' @param regime_ind A matrix of regime consistency indicators for a single
#'   regime.
#' @param regime A matrix of treatment assignments under the regime for each
#'   individual.
#' @param q_fits A list of fitted Q-function model objects for a single regime,
#'   one per stage.
#' @param k An integer indicating the stage index for which the derivative is
#'   computed.
#' @param feasible_sets_indicator A logical value indicating whether feasible sets
#'   are present in the trial design.
#' @param q_list A list of outcome regression model specifications (one per
#'   stage), used to determine the number of Q-function models at each stage
#'   (e.g., separate models for responders and non-responders).
#'
#' @return A numeric vector of partial derivatives of the value estimating
#'   equation with respect to the outcome regression parameters at stage
#'   \code{k}.
#'
#' @noRd
ee_dpsiv_dbeta <- function(df, pis, nus, regime_ind, regime, q_fits, k, feasible_sets_indicator, q_list){
  K <- length(nus$nu) - 1
  
  Ck <- t(apply(regime_ind,1,cumprod))
  CR <- rowSums(Ck)  # consistency by coarsening of regime. 
  # if someone is coarsened due to CR, then their level should be odd
  R <- pmin(CR*2 + 1, df$kappa*2)  # get the minimum coarsening level for each individual
  # instead of infinity, we use levels K*2+1, (K+1)*2
  
  # here, Beta_k is only used in r = 2*k-1 and 2*k terms. 
  # it also only in the augmentation terms, get design matrix here
  # Also, the design matrix should have modified ak to match regime k  
  df_mod <- df
  df_mod[,paste0('a',k)] <- regime[,k]
  
  if (length(q_fits[[k]]) > 1){
    designMatrix_r0 <- model.matrix(q_fits[[k]][["r0"]]@modelObj@model, data=df_mod)
    designMatrix_r1 <- model.matrix(q_fits[[k]][["r1"]]@modelObj@model, data=df_mod)
  } else {
    designMatrix <- model.matrix(q_fits[[k]]@modelObj@model, data=df_mod)
  }
  
  dPsi.dBetak <- lapply(vector(mode="list", length=nrow(df)), function(x) 0)
  
  for (r in (2*k-1):(2*k)){
    # now we get all the augmentation terms
    # the furthest stage reached with consistency is given as
    kr <- floor((r+1)/2) # here this should match k
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
    if (r %% 2 == 1){
      # odd r
      lambdar <- ( (pis[,kr]*(1-Ck[,kr]) + (1-pis[,kr])*Ck[,kr]) )  
      if (r==1){
        Kr <- ( (pis[,kr]*Ck[,kr] + (1-pis[,kr])*(1-Ck[,kr])) )  
      } else {
        Kr <- apply(pis[,1:(kr-1), drop = FALSE], 1, prod) * nus$nu[[kr]] * 
          ( (pis[,kr]*Ck[,kr] + (1-pis[,kr])*(1-Ck[,kr])) )  
      }
    } else{
      # even r
      # this is 0 for AIPW estimator 
      lambdar <-rep(1, nrow(df)) * (nus$nu[[kr]] - nus$nu[[kr+1]]) / nus$nu[[kr]]
      Kr <- apply(pis[,1:kr, drop = FALSE], 1, prod) * nus$nu[[kr+1]]
    }
    
    # change made 3/10/2022
    # the L functions (X^TBeta) may have Y for k=2 under feasible sets framework
    # this means that the derivative should be 0 if Y is used (i.e. R2=1 and kappai > 2)
    # the derivative should be non-zero at k=2 if we predicted X^tBeta using modified regime
    # this is determined by "new" indicators in q_step.R
    
    
    if (!feasible_sets_indicator){
      # multiple options for each stage every time
      new <- (df['kappa']>=k) # those who can attain a prediction from model  
      dOneTerm.dbetak <- lapply(X=1L:nrow(df), FUN = function(i) new[i] *
                                  ((R[i] >= r) * ((R[i] == r) - lambdar[i]) / Kr[i]) * designMatrix[i,])
      
    } else{# begin when feasible sets occurs
      if (!(paste0("r", k) %in% colnames(df))){ # are we at a stage where
        # there are Rk columns? If no, (this part), then treat as normal
        # multiple options for each stage every time
        new <- (df['kappa']>=k) # those who can attain a prediction from model  
        dOneTerm.dbetak <- lapply(X=1L:nrow(df), FUN = function(i) new[i] *
                                    ((R[i] >= r) * ((R[i] == r) - lambdar[i]) / Kr[i]) * designMatrix[i,])
        
      } else{ # there are rk at this stage
        # if this is true, then we need to determine if there are multiple 
        # models to fit. This occurs for interim analyses when we have some
        # individuals who have made it to stage k, had response determined,
        # but have not finished. 
        
        if (length(q_list[[k]]) > 1){
          ### start with those for rk==0 # non-responders ###
          new <- (df['kappa']>=k)*(df[paste0("r", k)] != 1) == 1
          dOneTerm.dbetak_r0 <- lapply(X=1L:nrow(df), FUN = function(i) new[i] *
                                      ((R[i] >= r) * ((R[i] == r) - lambdar[i]) / Kr[i]) * designMatrix_r0[i,])
          
          
          ### now do responders rk == 1 ###
          # we fit on k+1 individuals, but only use the model hats on stage k
          new <- (df['kappa']==k)*(df[paste0("r", k)] == 1) == 1
          dOneTerm.dbetak_r1 <- lapply(X=1L:nrow(df), FUN = function(i) new[i] *
                                      ((R[i] >= r) * ((R[i] == r) - lambdar[i]) / Kr[i]) * designMatrix_r1[i,])
          # we need to concatenate these since we want dPsiV/dBeta from 2 regressions at once
          # order R0 and R1 matches the ee_beta order
          dOneTerm.dbetak <- lapply(X=1L:nrow(df), FUN = function(i) c(dOneTerm.dbetak_r0[[i]], dOneTerm.dbetak_r1[[i]]))
          
        } else{ # this section may occur if we are estimating AIPW, no interim
          ### start with those for rk==0 # non-responders ###
          new <- (df['kappa']>=k)*(df[paste0("r", k)] != 1) == 1
          dOneTerm.dbetak <- lapply(X=1L:nrow(df), FUN = function(i) new[i] *
                                      ((R[i] >= r) * ((R[i] == r) - lambdar[i]) / Kr[i]) * designMatrix[i,])
          
        } # end handling q_list length;
      } # end the if/else for existence of responders w/i feasible sets
      
    } # end the if/else of feasible sets, 
    
    
    # dOneTerm.dbetak <- lapply(X=1L:nrow(df), FUN = function(i) 
    #   ((R[i] >= r) * ((R[i] == r) - lambdar[i]) / Kr[i]) * designMatrix[i,]) 
    
    # dOneTerm.dbetak <- lapply(X=1L:nrow(df), FUN = function(i) new[i] *
    #   ((R[i] >= r) * ((R[i] == r) - lambdar[i]) / Kr[i]) * designMatrix[i,])
     
    # end change 3/10/2022
    
    
    dPsi.dBetak <- lapply(X=1L:nrow(df), FUN = function(i) dPsi.dBetak[[i]] + dOneTerm.dbetak[[i]])
    
  } # end for loop  
  
  # we want to get dPsiV/dBetak
  return (Reduce("+", dPsi.dBetak))
}



#### With respect to V ####
# this should be +sum(Gamma_i) for all Vi, also is diagonal. 

#' Compute the derivative of the value estimating equation with respect to
#' the value parameters
#'
#' Computes \eqn{-\partial \Psi_V / \partial V}, which is a diagonal matrix
#' with entries equal to \eqn{-n_s} for each regime. This forms the lower-right
#' block of the An matrix.
#'
#' @param nregimes An integer indicating the number of treatment regimes being
#'   estimated.
#' @param nus A list as returned by \code{\link{get_nu}}, containing the
#'   estimated stage probabilities and the sample size \code{ns}.
#'
#' @return A diagonal matrix of dimension \code{nregimes} by \code{nregimes}
#'   with diagonal entries equal to \eqn{-n_s}.
#'
#' @noRd
ee_dpsiv_dv <- function(nregimes, nus){
  return (- nus$ns * diag(nregimes) )
}



