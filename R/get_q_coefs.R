#' Extract coefficients from all fitted Q-function models
#'
#' Extracts and concatenates the regression coefficients from all fitted
#' Q-function models across all regimes and stages.
#'
#' @param q_all A list of fitted outcome regression objects for each regime,
#'   as returned by \code{\link{estimate_values}}.
#'
#' @return A numeric vector of all Q-function regression coefficients, ordered
#'   by regime and then by stage.
#'
#' @export
get_q_coefs <- function(q_all){
  coefs <- c()
  for (ell in 1:length(q_all)){
    for (k in 1:length(q_all[[ell]]$q_fits)){
      if (length(q_all[[ell]]$q_fits[[k]]) > 1){
        coefs <- c(coefs,
                   unlist(lapply(q_all[[ell]]$q_fits[[k]], function(x) x@fitObj$coefficients)) )
      } else{
        coefs <- c(coefs,
                   q_all[[ell]]$q_fits[[k]]@fitObj$coefficients)
      }
    }

  }
  return (coefs)
}
