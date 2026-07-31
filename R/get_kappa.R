#' Compute the kappa (stage reached) for each individual
#'
#' Determines the total number of stages each individual has reached by a
#' specified analysis time \code{t_s}, based on their arrival times
#' \code{t1, ..., t_{K+1}}.
#'
#' @param df A data frame containing columns \code{t1, t2, ..., t_{K+1}}
#'   representing the times at which each individual reaches each stage.
#' @param t_s A numeric value specifying the analysis time point.
#' @param K An integer specifying the number of treatment stages (decision
#'   points).
#'
#' @return An integer vector of length \code{nrow(df)} indicating the number of
#'   stages each individual has reached by time \code{t_s}. Values range from 0
#'   (not yet enrolled) to \eqn{K+1} (outcome observed).
#'
#' @export
get_kappa <- function(df, t_s, K){
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame.")
  }
  if (!is.numeric(t_s) || length(t_s) != 1 || t_s < 0) {
    stop("'t_s' must be a non-negative numeric value.")
  }
  if (!is.numeric(K) || length(K) != 1 || K != as.integer(K) || K < 1) {
    stop("'K' must be a positive integer.")
  }
  return(rowSums(df[paste0('t', 1:(K+1))] <= t_s ))
}