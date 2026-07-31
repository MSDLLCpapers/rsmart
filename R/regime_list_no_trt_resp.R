
#' Generate regime lists for a two-stage SMART where responders are not
#' re-randomized
#'
#' Constructs the treatment assignment matrices and consistency indicator
#' matrices for each embedded regime in a two-stage SMART. For stages where
#' response status determines treatment deterministically, the regime matrix
#' is updated to reflect the responder treatment.
#'
#' @param emb_regimes A list of numeric vectors, each of length K, specifying
#'   the treatment codes (0 or 1) for non-responders at each stage. For
#'   example, `list(c(0,0), c(0,1), c(1,0), c(1,1))` encodes four regimes.
#' @param dat A data frame containing the SMART data. Expected to contain
#'   columns `a1, ..., aK` (treatment assignments) and optionally `r1, ...,
#'   rK` (response indicators). Typically this is the output from
#'   [gen_no_trt_resp()], but any data frame with the required columns may be
#'   used.
#' @param resp_trt A named list where each element corresponds to a stage with
#'   deterministic responder treatment (e.g., `"r2"`). Each element is a list
#'   of length equal to the number of regimes, specifying the treatment
#'   assignment (0 or 1) for responders at that stage under each regime.
#'   Default is `list("r2" = list(0, 0, 0, 0))`.
#'
#' @return A list of length `length(emb_regimes)`. Each element is a list with
#'   two components:
#'   \describe{
#'     \item{regime}{An n x K matrix of treatment assignments each individual
#'       would receive if they followed that regime.}
#'     \item{regime_ind}{An n x K indicator matrix (0 or 1) of whether each
#'       individual's observed treatment was consistent with the regime at each
#'       stage.}
#'   }
#'
#' @export
regime_list_no_trt_resp <- function(emb_regimes, dat, resp_trt=list("r2" = list(0, 0, 0, 0))){
  regime_list <- list()
  nRegimes <- length(emb_regimes)
  n <- nrow(dat)
  K <- length(emb_regimes[[1]])

  r <- 1
  for (regimen in emb_regimes){
    regime <- matrix(NA, nrow=n, ncol=K)
    regime_ind <- matrix(NA, nrow=n, ncol=K)

    for (k in 1:K){
      regime[,k] <- regimen[k]
      # update to reflect if a response at stage K is given
      if ((paste0("r", k) %in% colnames(dat)) & (!is.null(resp_trt))) {
        regime[,k] <- regime[,k]*(1 - dat[, paste0('r',k)]) +
          resp_trt[[paste0('r',k)]][[r]] * dat[, paste0('r',k)]
      }
      regime_ind[,k] <- (dat[, paste0('a',k)] == regime[,k])*1
    }
    colnames(regime) <- paste0('regime', 1:K)
    colnames(regime_ind) <- paste0('a', 1:K, '_ind')

    regime_list[[r]] <- list('regime' = regime,
                             "regime_ind" = regime_ind)
    r <- r+1
  }

  return(regime_list)
}
