#' Fit propensity score models for all stages
#'
#' Fits propensity score models at each stage of the SMART, using only
#' individuals who have reached that stage. Returns estimated propensity scores
#' and fitted model objects for all stages.
#'
#' @param df A data frame containing the trial data, including treatment
#'   assignments (\code{a1, a2, ...}), covariates, and a \code{kappa} column.
#' @param p_list A list of \code{modelObj} objects specifying the propensity
#'   score model at each stage.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{ps}{A data frame of estimated propensity scores with columns
#'       \code{pi1, ..., piK}. Individuals who have not reached a stage are
#'       assigned a value of 99.}
#'     \item{p_fits}{A list of fitted \code{modelObj} objects, one per stage.}
#'   }
#'
#' @export
pi_fits <- function(df, p_list){
  K <- length(x = p_list)
  fits <- list()
  ps <- data.frame(matrix(ncol=K, nrow=nrow(df)))
  # we want the names of the columns to be pi1,...,piK.
  colnames(ps) <- paste0("pi", c(1:K))

  for (k in 1:K){
    # only estimate for individuals with observed ak
    ones <- (df['kappa']>=k)

    pks <- pstep(pmodel = p_list[[k]],
                data = df[ones,, drop = FALSE],
                response = df[ones, paste0('a',k), drop = FALSE],
                k = k)

    fits[[k]] <- pks$pfit
    ps[,k] <- 99  # give everyone something.
    # these numbers do not factor in because the estimates pi should always be multiplied by the
    # indicator for kappa to ensure that it is 0 for anyone who has not reached
    # stage k. This is factored into the estimating equations.
    ps[ones, k] <- pks$ps  # only update those with estimates
  }

  return (list('ps' = ps,
               'p_fits' = fits))
}
