# now we have a data frame and we want to get all the terms that are applicable for individuals 
# We consider the 2K+1 coarsening terms here to keep code similar to paper
# We will do this differently from value_d2 and value_d3

# assume we have df, consistency indicators, Ck, qk, and pis at this point
# we also assume that kappa >= 1
# and we assume that Ci=1 => Delta=1
# We assume that the arrival times and treatment assignments are independent

#' Compute the individual-level value terms for a single regime
#'
#' Computes the \eqn{2K+1} coarsening-level value terms for each individual
#' under a single treatment regime. These include augmentation terms for levels
#' \eqn{r = 1, \ldots, 2K} and the IPW term for \eqn{R = \infty} (i.e.,
#' \eqn{R = 2K+1}). When outcome regression models (\code{qs}) are provided,
#' the augmented terms are computed; otherwise only the IPW term is non-zero.
#'
#' @param df A data frame containing the trial data, including treatment
#'   assignments, covariates, outcomes (\code{y}), and a \code{kappa} column.
#' @param regime_ind A matrix of regime consistency indicators for a single
#'   regime, with columns indicating whether each individual followed the
#'   regime at each stage.
#' @param pis A data frame of estimated propensity scores with columns
#'   \code{pi1, ..., piK}.
#' @param qs A list containing outcome regression fitted values for a single
#'   regime (as returned by \code{\link{get_q_fits}}), or \code{NULL} for IPW
#'   estimation.
#' @param nus A list as returned by \code{\link{get_nu}}, containing the
#'   estimated stage probabilities.
#'
#' @return A numeric matrix with \code{nrow(df)} rows and \eqn{2K+1} columns,
#'   where each column corresponds to a coarsening-level term in the value
#'   estimator.
#'
#' @export
value_terms <- function(df, regime_ind, pis, qs, nus){
  K <- length(nus$nu) - 1
    
  Ck <- t(apply(regime_ind,1,cumprod))
  CR <- rowSums(Ck)  # consistency by coarsening of regime. 
  # if someone is coarsened due to CR, then their level should be odd
  R <- pmin(CR*2 + 1, df$kappa*2)  # get the minimum coarsening level for each individual
  # instead of infinity, we use levels K*2+1, (K+1)*2
  
  # now we can get the augmentation terms for the correct levels of kappa
  
  # ipw term for R=infty
  Vterms <- matrix(0, nrow=nrow(df), ncol=2*K+1)
  rinfty <- df$y * (R == (2*K+1)) / (apply(pis, 1, prod) * nus$nu[[K+1]])
  Vterms[,2*K+1] <- rinfty
  
  if (!is.null(qs)){
    
    for (r in 1:(2*K)){
      # now we get all the augmentation terms
      # the furthest stage reached with consistency is given as
      k <- floor((r+1)/2)
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
      
      
      if (r %% 2 == 1){
        # odd r
        lambdar <- ( (pis[,k]*(1-Ck[,k]) + (1-pis[,k])*Ck[,k]) )  
        if (r==1){
          Kr <- ( (pis[,k]*Ck[,k] + (1-pis[,k])*(1-Ck[,k])) )  
        } else {
          Kr <- apply(pis[,1:(k-1), drop = FALSE], 1, prod) * nus$nu[[k]] * 
            ( (pis[,k]*Ck[,k] + (1-pis[,k])*(1-Ck[,k])) )  
        }
      } else{
        # even r
        # this is 0 for AIPW estimator 
        lambdar <- (nus$nu[[k]] - nus$nu[[k+1]]) / nus$nu[[k]]
        Kr <- apply(pis[,1:k, drop = FALSE], 1, prod) * nus$nu[[k+1]]
      }
      
      # here the Q should be using the modified regimes
      # first, make sure to only take the terms for which R >=r
      augterm <- ifelse((R >= r), (((R == r) - lambdar) / Kr) * qs$mod_regime_vhats[,k], 0)
      # next, correct for when Kr = 0 due to regime inconsistency at the current stage
      # augterm <- ifelse(Ck[,k] == 0, 0, augterm)
      
      Vterms[,r] <- augterm
      
      # update Kr lag
      Krlag1 <- Kr
    }
    
  }
  
  return(Vterms)
}
