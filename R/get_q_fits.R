#' Fit outcome regression (Q-function) models across all stages for a single
#' regime
#'
#' Fits Q-function models backwards from the last stage to the first for a
#' single treatment regime. At each stage, both regime-modified and unmodified
#' predicted values are computed. Handles feasible sets (where responders may
#' not be re-randomized) and interim analyses.
#'
#' @param df A data frame containing the trial data, including treatment
#'   assignments, covariates, outcomes (\code{y}), and a \code{kappa} column.
#' @param q_list A list of outcome regression model specifications (one per
#'   stage). Each element is either a \code{modelObj} object or a named list
#'   with elements \code{"r0"} and \code{"r1"} for non-responders and
#'   responders, respectively.
#' @param regime A matrix of treatment assignments under the regime being
#'   evaluated, with columns corresponding to stages.
#' @param feasible_sets_indicator A logical value indicating whether feasible sets
#'   are present. Default is \code{FALSE}. If individuals are not re-randomized
#'  (i.e. for some stage a response status prevents individuals from receiving a
#'  random treatment), then feasible_sets_indicator should be set to TRUE.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{q_fits}{A list of fitted \code{modelObj} objects, one per stage.
#'       When feasible sets with multiple models are used, the element is a
#'       named list with \code{"r0"} and \code{"r1"}.}
#'     \item{mod_regime_vhats}{A matrix of regime-modified predicted values
#'       with columns \code{q1, ..., q_{K+1}}.}
#'     \item{unmod_regime_vhats}{A matrix of unmodified predicted values with
#'       columns \code{q1_nochange, ..., q_{K+1}_nochange}.}
#'   }
#'
#' @export
get_q_fits <- function(df, q_list, regime, feasible_sets_indicator=FALSE){
  # total number of decision points
  K <- length(q_list)
  # create storage for modified and unmodified predictors
  unmod_regime_vhats <- matrix(data = 0.0, nrow = nrow(df), ncol = K + 1)
  mod_regime_vhats <- matrix(data = 0.0, nrow = nrow(df), ncol = K + 1)

  # for Kappa = decisions + 1 keep the response
  unmod_regime_vhats[,K + 1] <- df$y
  mod_regime_vhats[,K + 1] <- df$y

  q_fits <- list()

  for (k in (K:1)){
    # feasible sets indicator is a boolean that indicates if some stages have
    # only 1 treatmenat (TRUE) or multiple (FALSE)

    # pass previous expected outcome back
    # these will be updated by regimes in following code
    mod_regime_vhats[, k] <- mod_regime_vhats[,k+1]
    unmod_regime_vhats[, k] <- unmod_regime_vhats[,k+1]


    if (!feasible_sets_indicator){
      # multiple options for each stage every time
      ones <- (df['kappa']>=(k+1)) # those who have finished
      new <- (df['kappa']>=k) # those who can attain a prediction from model

      # now fit the model and get predictions
      vk <- qstep(qmodel = q_list[[k]],
                  data = df[ones,, drop = FALSE],
                  response = mod_regime_vhats[ones, k+1, drop = FALSE],
                  newdata = df[new,, drop = FALSE],
                  regime = regime[new,k],
                  txName = paste0('a',k) )

      # update for those who were able to receive
      mod_regime_vhats[new, k] <- vk$hats_mod
      unmod_regime_vhats[new, k] <- vk$hats_unmod

      # add the qfit to the outlist
      q_fits[[k]] <- vk$qfit

    } else{# begin when feasible sets occurs
      if (!(paste0("r", k) %in% colnames(df))){ # are we at a stage where
        # there are Rk columns? If no, (this part), then treat as normal
        # multiple options for each stage every time
        ones <- (df['kappa']>=(k+1)) # those who have finished
        new <- (df['kappa']>=k) # those who can attain a prediction from model

        # now fit the model and get predictions
        vk <- qstep(qmodel = q_list[[k]],
                    data = df[ones,, drop = FALSE],
                    response = mod_regime_vhats[ones, k+1, drop = FALSE],
                    newdata = df[new,, drop = FALSE],
                    regime = regime[new,k],
                    txName = paste0('a',k) )

        # update for those who were able to receive
        mod_regime_vhats[new, k] <- vk$hats_mod
        unmod_regime_vhats[new, k] <- vk$hats_unmod

        # add the qfit to the outlist
        q_fits[[k]] <- vk$qfit
      } else{ # there are rk at this stage
        # if this is true, then we need to determine if there are multiple
        # models to fit. This occurs for interim analyses when we have some
        # individuals who have made it to stage k, had response determined,
        # but have not finished.


        if (length(q_list[[k]]) > 1){
          # this means that we have multiple functions to fit
          # q_fits[[k]] <- list("r0", "r1")

          ### start with those for rk==0 # non-responders ###
          ones <- (df['kappa']>=(k+1))*(df[paste0("r", k)] != 1) == 1
          new <- (df['kappa']>=k)*(df[paste0("r", k)] != 1) == 1

          # now fit the model and get predictions
          # note that "response" in qstep is the outcome, not trt response
          vk <- qstep(qmodel = q_list[[k]][[paste0("r", 0)]],
                      data = df[ones,, drop = FALSE],
                      response = mod_regime_vhats[ones, k+1, drop = FALSE],
                      newdata = df[new,, drop = FALSE],
                      regime = regime[new,k],
                      txName = paste0('a',k) )

          # update for those who were able to receive
          mod_regime_vhats[new, k] <- vk$hats_mod
          unmod_regime_vhats[new, k] <- vk$hats_unmod

          # add the qfit to the outlist
          qr0 <- vk$qfit


          ### now do responders rk == 1 ###
          # we fit on k+1 individuals, but only use the model hats on stage k
          ones <- (df['kappa']==(k+1))*(df[paste0("r", k)] == 1) == 1
          new <- (df['kappa']==k)*(df[paste0("r", k)] == 1) == 1

          # now fit the model and get predictions
          # note that "response" in qstep is the outcome, not trt response
          vk <- qstep(qmodel = q_list[[k]][[paste0("r", 1)]],
                      data = df[ones,, drop = FALSE],
                      response = mod_regime_vhats[ones, k+1, drop = FALSE],
                      newdata = df[new,, drop = FALSE],
                      regime = regime[new,k],
                      txName = paste0('a',k) )

          # update for those who were able to receive
          mod_regime_vhats[new, k] <- vk$hats_mod
          unmod_regime_vhats[new, k] <- vk$hats_unmod

          # add the qfit to the outlist
          qr1 <- vk$qfit

          q_fits[[k]] <- list("r0" = qr0,
                              "r1" = qr1)

        } else{ # this section may occur if we are estimating AIPW, no interim
          ### start with those for rk==0 # non-responders ###
          ones <- (df['kappa']>=(k+1))*(df[paste0("r", k)] != 1) == 1
          new <- (df['kappa']>=k)*(df[paste0("r", k)] != 1) == 1

          # now fit the model and get predictions
          # note that "response" in qstep is the outcome, not trt response
          vk <- qstep(qmodel = q_list[[k]],
                      data = df[ones,, drop = FALSE],
                      response = mod_regime_vhats[ones, k+1, drop = FALSE],
                      newdata = df[new,, drop = FALSE],
                      regime = regime[new,k],
                      txName = paste0('a',k) )

          # update for those who were able to receive
          mod_regime_vhats[new, k] <- vk$hats_mod
          unmod_regime_vhats[new, k] <- vk$hats_unmod

          # add the qfit to the outlist
          q_fits[[k]] <- vk$qfit

        } # end handling q_list length;
      } # end the if/else for existence of responders w/i feasible sets


    } # end the if/else of feasible sets,

  } # end for loop

  colnames(mod_regime_vhats) <- paste0('q', 1:(K+1))
  colnames(unmod_regime_vhats) <- paste0(paste0('q', 1:(K+1)), '_nochange')
  return(list('q_fits' = q_fits,
              'mod_regime_vhats' = mod_regime_vhats,
              'unmod_regime_vhats' = unmod_regime_vhats))
} # end get_q_fits
