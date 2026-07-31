#' Compute individual-level estimating equations for propensity score parameters
#' (pi)
#'
#' Computes \eqn{\Psi_{\pi,i}} for each individual \eqn{i}, the estimating
#' equation contributions for the propensity score model parameters across all
#' stages. These are used as part of the Bn matrix in the sandwich variance
#' estimator.
#'
#' @param df A data frame containing the trial data, including treatment
#'   assignments (\code{a1, a2, ...}), covariates, and a \code{kappa} column.
#' @param p_fits A list of fitted propensity score model objects (one per stage),
#'   each a \code{modelObj} fit object.
#'
#' @return A numeric matrix with rows corresponding to individuals and columns
#'   corresponding to propensity score parameters across all stages, representing
#'   the individual-level estimating equation contributions for pi.
#'
#' @noRd
ee_psi_pi <- function(df, p_fits){
  g <- function(xi, thetapi){
    # xi is 1x4, thetapi is a vector of length 4
    return (exp(xi %*% thetapi) / (1 + exp(xi %*% thetapi)))
  }
  gprime <- function(xi, thetapi){
    # note that xi is 1x4 and the expfn is 1x1, result should be 4x1
    return (t(xi) %*% (exp(xi %*% thetapi) / ((1 + exp(xi %*% thetapi))^2)))
  }
  # second derivative of g
  gprime2 <- function(xi, thetapi){
    # xi is 1x4
    return (t(xi) %*% ( exp(xi %*% thetapi) / ((1 + exp(xi %*% thetapi))^3) ) %*% xi)
  }
  # used to format the output from lapply correctly
  # using t(sapply(...)) does not work if there is only 1 model parameter
  makeNice <- function(...) {matrix(unlist(...), nrow=length(...), byrow=TRUE)}
  
  
  psi.pi <- c()
  for (k in c(1:length(p_fits))){
    
    k_ind <- (df[,'kappa']>=k)*1
    dfkm <- model.matrix(p_fits[[k]]@modelObj@model, data=df)
    
    psii <- makeNice(
      lapply(X=1L:nrow(df), FUN = function(i) drop(df[,paste0('a',k)][i] - 
                                                     g(dfkm[i,,drop = FALSE], 
                                                       p_fits[[k]]@fitObj$coefficients)) *
               k_ind[i] *
               gprime(dfkm[i,,drop = FALSE], p_fits[[k]]@fitObj$coefficients)
      )
    )
    
    
    psi.pi <- cbind(psi.pi, psii) 
  }
  return(psi.pi)
}


#' Compute the negative derivative of the propensity score estimating equations
#' with respect to pi
#'
#' Computes \eqn{-\partial \Psi_{\pi} / \partial \pi}, the negative derivative
#' of the propensity score estimating equations with respect to the propensity
#' score model parameters. This forms a block diagonal component of the An
#' matrix in the sandwich variance estimator. Based on Equation 7.29 in Boos
#' and Stefanski (p. 319).
#'
#' @param df A data frame containing the trial data, including treatment
#'   assignments (\code{a1, a2, ...}), covariates, and a \code{kappa} column.
#' @param p_fits A list of fitted propensity score model objects (one per stage),
#'   each a \code{modelObj} fit object.
#'
#' @return A block diagonal matrix (potentially sparse via \code{Matrix::bdiag})
#'   representing \eqn{-\partial \Psi_{\pi} / \partial \pi} across all stages.
#'
#' @noRd
ee_dpsi_pi <- function(df, p_fits){
  g <- function(xi, thetapi){
    # xi is 1x4, thetapi is a vector of length 4
    return (exp(xi %*% thetapi) / (1 + exp(xi %*% thetapi)))
  }
  gprime <- function(xi, thetapi){
    # note that xi is 1x4 and the expfn is 1x1, result should be 4x1
    return (t(xi) %*% (exp(xi %*% thetapi) / ((1 + exp(xi %*% thetapi))^2)))
  }
  # second derivative of g
  gprime2 <- function(xi, thetapi){
    # xi is 1x4
    return (t(xi) %*% ( exp(xi %*% thetapi) / ((1 + exp(xi %*% thetapi))^3) ) %*% xi)
  }
  # used to format the output from lapply correctly
  # using t(sapply(...)) does not work if there is only 1 model parameter
  # makeNice <- function(...) {matrix(unlist(...), nrow=length(...), byrow=TRUE)}
  
  # general formula is given as eq 7.29 pg 319 Boos Stefanski
  # \sum_{i=1}^N gprime gprime^t - {(Y_i - g)g''}
  for (k in c(1:length(p_fits))){
    
    k_ind <- (df[,'kappa']>=k)*1
    dfkm <- model.matrix(p_fits[[k]]@modelObj@model, data=df)
    
    dpsi.pi_k <- Reduce("+", 
                        lapply(X=1L:nrow(df), FUN = function(i) k_ind[i] * 
                                 gprime(dfkm[i,,drop = FALSE], 
                                        p_fits[[k]]@fitObj$coefficients) %*% 
                                 t(gprime(dfkm[i,,drop = FALSE], 
                                          p_fits[[k]]@fitObj$coefficients)) 
                               # remove because expectation is mean 0, and adds computational burden. 
                               # - (df[i,'kappa']>=k) * {
                               #   drop(df[,paste0('a',k)][i] - 
                               #          g(model.matrix(p_fits[[k]]@modelObj@model, data=df[i,]), 
                               #            p_fits[[k]]@fitObj$coefficients)) * 
                               #     gprime2(model.matrix(p_fits[[k]]@modelObj@model, data=df[i,]), p_fits[[k]]@fitObj$coefficients)
                               # }
                        )
    )
    # this should be block diagonal
    if (k == 1){
      dpsi.pi <- dpsi.pi_k
    } else{
      dpsi.pi <- Matrix::bdiag(dpsi.pi, dpsi.pi_k)
    }
  }
  return(-dpsi.pi)
  return(-dpsi.pi)
}


    

