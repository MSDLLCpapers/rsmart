#' Estimate values for all treatment regimes
#'
#' Wraps the outcome regression and value term functions to estimate the value
#' of all regimes in the provided list. For each regime, Q-functions are fitted
#' (if specified), augmentation and IPW terms are computed, and the regime value
#' is estimated.
#'
#' @param df A data frame containing the trial data, including treatment
#'   assignments, covariates, outcomes, and a \code{kappa} column.
#' @param q_list A list of outcome regression model specifications (one per
#'   stage), or \code{NULL} for IPW estimation.
#' @param regime_all A list of regime objects, each containing a \code{regime}
#'   matrix and a \code{regime_ind} indicator matrix.
#' @param feasible_sets_indicator A logical value indicating whether feasible sets
#'   are present in the trial design.
#' @param pis A data frame of estimated propensity scores with columns
#'   \code{pi1, ..., piK}.
#' @param nus A list as returned by \code{\link{get_nu}}, containing the
#'   estimated stage probabilities.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{value}{A numeric vector of estimated regime values.}
#'     \item{df}{A list of matrices of value term components, augmentation terms 1
#'       through 2K then the IPW term for each regime.}
#'     \item{q_all}{A list of fitted Q-function objects for each regime, or an
#'       empty list if \code{q_list} is \code{NULL}.}
#'   }
#'
#' @export
estimate_values <- function(df, q_list, regime_all, feasible_sets_indicator, 
                            pis, nus){
  # get value estimate
  values <- c()
  vTermDfs <- list()
  q_all <- list()
  L <- length(regime_all)
  for (ell in 1:L){
    if (!is.null(q_list)){
      qs <- get_q_fits(df, q_list, regime_all[[ell]]$regime, feasible_sets_indicator=feasible_sets_indicator)  
    } else{
      qs <- NULL
    }
    vTerms <- value_terms(df, regime_all[[ell]]$regime_ind, pis, qs, nus)
    value <- sum(vTerms) / sum(df$kappa > 0)
    vTermDfs[[ell]] <- vTerms
    values <- c(values, value)
    q_all[[ell]] <- qs  # if qs is null this will not assign anything
  }
  return(list('value' = values,
              'df' = vTermDfs,
              'q_all' = q_all))
}


