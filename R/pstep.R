#' Fit a propensity score model for a single stage
#'
#' Fits a propensity score model at a single stage of the SMART and returns the
#' fitted probabilities. The probability that each individual received the
#' treatment they actually received is computed.
#'
#' @param pmodel A \code{modelObj} object specifying the propensity score model
#'   (typically a binomial GLM).
#' @param data A data frame of individuals to fit the model on (typically those
#'   who have reached stage \code{k}).
#' @param response A vector or data frame column of treatment assignments at
#'   stage \code{k}.
#' @param k An integer indicating the stage number.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{pk}{A numeric vector of predicted probabilities that treatment
#'       is 1 at stage \code{k}.}
#'     \item{ps}{A numeric vector of estimated propensity scores (probability
#'       that the individual received the treatment they actually received).}
#'     \item{pfit}{The fitted \code{modelObj} object.}
#'   }
#'
#' @export
pstep <- function(pmodel, data, response, k) {
  # fit the qstep model
  pfit <- modelObj::fit(object = pmodel, data = data, response = response)
  # get unmodified predicted values for newdata
  pk <- modelObj::predict(object = pfit, newdata = data)
  # these are now the probability that the outcome is 1
  # however, we want the probability that individual receives the
  # treatment that they did indeed get
  ps <- {pk*{ data[,paste0('a',k)] == 1L} + {1.0 - pk}*{data[,paste0('a',k)] == 0L}}
  return (list('pk' = pk,
               'ps' = ps,
               'pfit' = pfit))
}
