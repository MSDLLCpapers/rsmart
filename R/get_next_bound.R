#' Compute subsequent stopping boundaries for a group sequential design
#'
#' Given the boundary from the first analysis, determines the stopping boundary
#' for the \eqn{s}-th analysis that controls the overall type I error rate.
#' Uses an iterative search with the joint distribution of test statistics
#' across analyses.
#'
#' @param alpha A numeric value specifying the overall type I error rate to
#'   control. Default is 0.05.
#' @param inf_frac A numeric vector of information fractions indicating when
#'   analyses are conducted.
#' @param spend_fn A character string specifying the alpha spending function.
#'   Either \code{"OF"} (O'Brien-Fleming) or \code{"Pocock"}.
#' @param corr A correlation matrix of Z-statistics across all analyses and
#'   regimes. Should have dimension \eqn{SL \times SL} by \eqn{SL \times SL}.
#' @param test_type A character string specifying the type of test to be
#'   performed. Either \code{"one-sided"} or \code{"two-sided"}. Default is
#'   \code{"one-sided"}.
#' @param lambda A numeric value for the initial step size used in the
#'   iterative boundary search. Default is 0.1.
#' @param tol A numeric value specifying the convergence tolerance. Default is
#'   1e-6.
#' @param prev_bound A numeric vector of previously computed boundaries from
#'   analyses \eqn{1, \ldots, s-1}, with dimension \eqn{(s-1) \times L}.
#' @param s An integer indicating the current analysis number. Default is 2.
#' @param max_iter A positive integer specifying the maximum number of
#'   iterations. Default is 1000.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{bound}{The computed stopping boundary for analysis \code{s}.}
#'     \item{typeI}{The achieved cumulative type I error rate through analysis
#'       \code{s}.}
#'     \item{convergence}{A logical value indicating whether the algorithm
#'       converged.}
#'     \item{iters}{The number of iterations used.}
#'   }
#'
#' @export
get_next_bound <- function(alpha = 0.05, inf_frac = c(0.5, 1), spend_fn = "OF",
                          corr = diag(x=1, nrow=1, ncol=1), test_type = "one-sided",
                          lambda = 0.1, tol = 1e-6, prev_bound = NULL, s = 2,
                         max_iter = 1000){
  # previous bound should be dimension (s-1)*L
  
  S <- length(inf_frac)  # number of analyses to perform
  L <- dim(corr)[1] / S  # get number of regimes estimated each time point 
  
  inf <- inf_frac[s]
  # begin by finding the first boundary 
  # note that the superiority test uses 1-alpha. 
  # to change this for two-sided testing, use 1-(alpha/2)
  if (spend_fn == "OF"){
    alphat <- 2 - 2 * pnorm(qnorm(p=1-alpha/2, lower.tail=TRUE) / sqrt(inf) )
  } else if (spend_fn == "Pocock") {
    alphat <- alpha * log(1 + (exp(1) -1 ) * inf)
  }
  
  corrs <- corr[1:(s*L),
                1:(s*L)]

  # we want to choose a starting boundary that will inflate the type I error rate so we can use the same code as before...
  # if we select the boundary as the boundary for a single analysis single regime type I error rate, it will be too low, 
  # so we start with qnorm(1-alpha)
  bs <- qnorm(1-alpha)
  
  bounds <- c(prev_bound, rep(bs, L))
  
  if (test_type == "two-sided") {
    lower = c(-prev_bound, rep(-bs, L))
  } else {
    lower = -Inf
  }

  typeI_ <- 1- pmvnorm(lower=lower, 
                       upper=bounds, 
                       mean=rep(0, nrow(corrs)), 
                       corr = corrs )
  typeI_  # we anticipate that the type I error rate is inflated, so b1 must increase 
  
  alphatarget <- alpha 
  
  loopcount <- 0
  converged <- TRUE
  
  while(TRUE){  # loop to update lambda by a factor of 10
    while (TRUE) {  # loop to change bound by lambda until over alpha
      loopcount <- loopcount + 1
        if (test_type == "two-sided") {
          lower = c(-prev_bound, rep(-bs, L))
        }
      typeI <- 1- pmvnorm(lower=lower, 
                          upper= c(prev_bound, rep(bs, L)),
                          mean=rep(0, nrow(corrs)),
                          corr = corrs )
      if (typeI < alphatarget){
        break
      }
      # update the boundaries based on which side of alpha we are on
      else{
        # update the boundaries
        bs <- bs + lambda
      }
      
      if (loopcount > max_iter){
        converged <- FALSE
        break
      }
    }
    # break if within tolerable window
    if ((typeI <= alphatarget + tol) & (typeI >= alphatarget - tol)){
      break
    }
    # else update the lambda 
    else{
      bs <- bs - lambda
      lambda <- lambda / 10 
    }
    
    if (loopcount > max_iter){
      converged <- FALSE
      break
    }
  }
  return(list("bound" = bs,
              "typeI" = typeI,
              "convergence" = converged,
              "iters" = loopcount))
    
}
