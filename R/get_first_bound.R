#' Compute the first stopping boundary for a group sequential design
#'
#' Determines the first analysis stopping boundary that controls the familywise
#' type I error rate at a specified level, adjusting for the multiplicity of
#' multiple treatment regimes. Uses an iterative search to find the boundary
#' value. This function works only for superiority (one-sided) testing.
#'
#' @param alpha A numeric value specifying the overall type I error rate to
#'   control. Default is 0.05.
#' @param inf_frac A numeric vector of information fractions indicating when
#'   analyses are conducted. Values should be between 0 and 1.
#' @param spend_fn A character string specifying the alpha spending function.
#'   Either \code{"OF"} (O'Brien-Fleming) or \code{"Pocock"}.
#' @param corr A correlation matrix of Z-statistics at the first time point,
#'   accounting for multiple regimes. Default is a 1x1 identity matrix.
#' @param test_type A character string specifying the type of test to be
#'   performed. Either \code{"one-sided"} or \code{"two-sided"}. Default is
#'   \code{"one-sided"}.
#' @param lambda A numeric value for the initial step size used in the
#'   iterative boundary search. Default is 0.1.
#' @param tol A numeric value specifying the convergence tolerance for the
#'   type I error boundary. Default is 1e-6.
#' @param max_iter A positive integer specifying the maximum number of
#'   iterations for the boundary search. Default is 1000.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{bound}{The computed stopping boundary for the first analysis.}
#'     \item{typeI}{The achieved type I error rate at the boundary.}
#'     \item{convergence}{A logical value indicating whether the algorithm
#'       converged.}
#'     \item{iters}{The number of iterations used.}
#'   }
#'
#' @export
get_first_bound <- function(alpha = 0.05, inf_frac = c(0.5, 1), spend_fn = "OF",
                          corr = diag(x=1, nrow=1, ncol=1), test_type = "one-sided",
                          lambda = 0.1, tol = 1e-6, max_iter = 1000){
  inf <- inf_frac[1]
  # begin by finding the first boundary 
  # alphat is the total alpha spent at analysis with inf

  # note that the superiority test uses 1-alpha. 
  if (spend_fn == "OF"){
    alphat <- 2 - 2 * pnorm(qnorm(p=1-alpha/2, lower.tail=TRUE) / sqrt(inf) )
  } else if (spend_fn == "Pocock") {
    alphat <- alpha * log(1 + (exp(1) -1 ) * inf )
  }
  # now we adjust the boundary for the multiplicity of multiple regimes 
  # so that Pr(Z(t_1) > c_1) = \alpha(t_1)
  # the initial boundary will be too small so we have to increase the bound
  if (test_type == "two-sided") {
    c1 <- qnorm(1-alphat/2)
    b1 <- c1
    lower = rep(-b1, nrow(corr))
  } else{
    lower = -Inf
    c1 <- qnorm(1-alphat)
    b1 <- c1
  }
  typeI_ <- 1- pmvnorm(lower=lower, 
                       upper=rep(b1, nrow(corr) ), 
                       mean=rep(0, nrow(corr)), 
                       corr = corr)

  alphatarget <- alphat 
  
  loopcount <- 0
  converged <- TRUE
  
  while(TRUE){  # loop to update lambda by a factor of 10
    while (TRUE) {  # loop to change bound by lambda until over alpha
      loopcount <- loopcount + 1
        if (test_type == "two-sided") {
          lower = rep(-b1, nrow(corr))
        }
        typeI <- 1 - pmvnorm(lower=lower, 
                          upper=rep(b1, nrow(corr) ), 
                          mean=rep(0, nrow(corr)), 
                          corr = corr )
      if (typeI < alphatarget){
        break
      }
      # update the boundaries based on which side of alpha we are on
      else{
        # update the boundaries
        b1 <- b1 + lambda
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
      b1 <- b1 - lambda
      lambda <- lambda / 10 
    }
    
    if (loopcount > max_iter){
      converged <- FALSE
      break
    }
  }
  # boundary is b1, with type I error typeI
  return(list("bound" = b1,
              "typeI" = typeI,
              "convergence" = converged,
              "iters" = loopcount))
}
