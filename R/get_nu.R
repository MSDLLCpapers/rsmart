#' Estimate stage arrival probabilities (nu)
#'
#' Estimates the probability that an individual has reached each stage given
#' they are enrolled in the trial. Requires the data frame to have a
#' \code{kappa} column indicating the stage reached. The returned list uses
#' indexing such that \code{nu[[k]]} corresponds to stage \eqn{k}, and
#' \code{nu[[K+1]]} is the probability that an individual has their final
#' outcome observed.
#'
#' @param df A data frame containing a \code{kappa} column indicating the
#'   stage reached by each individual.
#' @param K An integer specifying the number of treatment stages (decision
#'   points).
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{nu}{A list of length \eqn{K+1} where \code{nu[[k]]} is the
#'       estimated probability that an individual has reached stage \eqn{k}
#'       given enrollment.}
#'     \item{ns}{The total number of individuals enrolled in the trial
#'       (\code{kappa > 0}).}
#'     \item{nd}{The proportion of individuals who have reached their last
#'       treatment stage and have their outcome observed.}
#'   }
#'
#' @export
get_nu <- function(df, K){
  if (!("kappa" %in% colnames(df))) {
    stop("'df' must contain a 'kappa' column.")
  }
  ns <- sum((df['kappa'] > 0))
  if (ns <= 0) {
    stop("No enrolled individuals found (ns must be a positive integer).")
  }
  nu <- list()
  for (k in ((K+1):1)){
    nu[[k]] <- sum((df['kappa']>=k)*((df['kappa'] > 0))) / ns
  }
  # proportion of those who have received their last treatment,
  # but have not finished
  if (sum(df['kappa'] >= K) == 0) {
    stop("No individuals reached stage K (nd denominator is zero).")
  }
  nd <- sum(df['kappa']>=(K+1)) / sum(df['kappa'] >= K)
  return(list('nu' = nu,
              'ns' = ns,
              'nd' = nd))
}
