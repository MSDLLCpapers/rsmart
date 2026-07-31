#' Fit an outcome regression (Q-function) model for a single stage
#'
#' Fits a single Q-function model at one stage of the SMART. Predictions are
#' obtained both for the unmodified data and for data where the treatment is set
#' to the recommended regime.
#'
#' @param qmodel A \code{modelObj} object specifying the outcome regression
#'   model.
#' @param data A data frame of individuals used to fit the model (typically
#'   those who have completed the stage).
#' @param response A numeric vector or single-column data frame of responses
#'   (pseudo-outcomes from later stages or observed outcomes).
#' @param newdata A data frame of individuals for whom predictions are desired
#'   (typically those who have reached the stage).
#' @param regime A vector of recommended treatment assignments under the regime
#'   being evaluated.
#' @param txName A character string specifying the name of the treatment column
#'   to modify (e.g., \code{"a1"}, \code{"a2"}).
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{hats_mod}{A numeric vector of predicted values under the
#'       regime-consistent treatment.}
#'     \item{hats_unmod}{A numeric vector of predicted values under the
#'       actual (unmodified) treatment.}
#'     \item{qfit}{The fitted \code{modelObj} object.}
#'   }
#'
#' @export
qstep <- function (qmodel, data, response, newdata, regime, txName) {
  # fit the qstep model
  qfit <- modelObj::fit(object = qmodel, data = data, response = response)
  # get unmodified predicted values for newdata
  hats_unmod <- modelObj::predict(object = qfit, newdata = newdata)
  # set tx to recommended
  newdata[,txName] <- regime
  # predict the newdata with regime consistent estimates
  hats_mod <- modelObj::predict(object = qfit, newdata = newdata)
  return (list('hats_mod' = hats_mod,
               'hats_unmod' = hats_unmod,
               'qfit' = qfit))
}
