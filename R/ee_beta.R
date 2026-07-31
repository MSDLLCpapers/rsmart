# Functions to support computing the estimating equations for Beta
# to be used in sandwich variance estimation corresponding to the q functions.
# We expect qs to be the list of q functions for a single regime.
# We perform this for all ell estimated regimes.

#' Compute individual-level estimating equations for outcome regression
#' coefficients (beta)
#'
#' Computes \eqn{\Psi_{\beta,i}} for each individual \eqn{i}, the estimating
#' equation contributions for the outcome regression (Q-function) parameters
#' across all regimes and stages. These are used as part of the Bn matrix in the
#' sandwich variance estimator.
#'
#' @param df A data frame containing the trial data, including treatment
#'   assignments, covariates, outcomes, and a \code{kappa} column.
#' @param q_all A list of fitted outcome regression objects for each regime,
#'   as returned by \code{\link{estimate_values}}.
#' @param feasible_sets_indicator A logical value indicating whether feasible sets
#'   are present in the trial design.
#' @param q_list A list of outcome regression model specifications (one per
#'   stage), used to determine the number of Q-function models at each stage
#'   (e.g., separate models for responders and non-responders).
#'
#' @return A numeric matrix with rows corresponding to individuals and columns
#'   corresponding to outcome regression parameters across all regimes and
#'   stages, representing the individual-level estimating equation contributions
#'   for beta.
#'
#' @noRd
ee_psi_beta <- function(df, q_all, feasible_sets_indicator, q_list){
  psi.beta <- c()
  for (ell in c(1:length(q_all))){
    psi.beta.ell <- c()
    for (k in c(1:length(q_all[[ell]]$q_fits))){
      # get the model matrix
      # first check that there are feasible sets at time step k
      if (!feasible_sets_indicator){
        # multiple options for each stage every time
        # no feasible sets
        ones <- (df['kappa']>=(k+1)) # those who have finished
        
        qk_design <- model.matrix(q_all[[ell]]$q_fits[[k]]@modelObj@model, data=df)  
        
        psi.betak <- diag(c(q_all[[ell]]$mod_regime_vhats[,k+1]*ones - q_all[[ell]]$unmod_regime_vhats[,k]*ones)) %*%
          qk_design
        psi.beta.ell <- cbind(psi.beta.ell, psi.betak) 

      } else{# begin when feasible sets occurs
        if (!(paste0("r", k) %in% colnames(df))){ # are we at a stage where
          # there are Rk columns? If no, (this part), then treat as normal
          # multiple options for each stage every time
          ones <- (df['kappa']>=(k+1)) # those who have finished
          
          qk_design <- model.matrix(q_all[[ell]]$q_fits[[k]]@modelObj@model, data=df)  
          
          psi.betak <- diag(c(q_all[[ell]]$mod_regime_vhats[,k+1]*ones - q_all[[ell]]$unmod_regime_vhats[,k]*ones)) %*%
            qk_design
          psi.beta.ell <- cbind(psi.beta.ell, psi.betak) 
          
          
        } else{ # there are rk at this stage
          # if this is true, then we need to determine if there are multiple 
          # models to fit. This occurs for interim analyses when we have some
          # individuals who have made it to stage k, had response determined,
          # but have not finished. 
          
          if (length(q_list[[k]]) > 1){
            # this means that we have multiple functions to fit
            
            ### start with those for rk==0 # non-responders ###
            ones <- (df['kappa']>=(k+1))*(df[paste0("r", k)] != 1) == 1
            
            qk_design <- model.matrix(q_all[[ell]]$q_fits[[k]][["r0"]]@modelObj@model, data=df)  
            
            psi.betak <- diag(c(q_all[[ell]]$mod_regime_vhats[,k+1]*ones - q_all[[ell]]$unmod_regime_vhats[,k]*ones)) %*%
              qk_design
            psi.beta.ell <- cbind(psi.beta.ell, psi.betak) 
            
            ### now do responders rk == 1 ###
            # we fit on k+1 individuals, but only use the model hats on stage k
            ones <- (df['kappa']==(k+1))*(df[paste0("r", k)] == 1) == 1
            
            qk_design <- model.matrix(q_all[[ell]]$q_fits[[k]][["r1"]]@modelObj@model, data=df)  
            
            psi.betak <- diag(c(q_all[[ell]]$mod_regime_vhats[,k+1]*ones - q_all[[ell]]$unmod_regime_vhats[,k]*ones)) %*%
              qk_design
            psi.beta.ell <- cbind(psi.beta.ell, psi.betak) 
            
          } else{ # this section may occur if we are estimating AIPW, no interim
            ### start with those for rk==0 # non-responders ###
            ones <- (df['kappa']>=(k+1))*(df[paste0("r", k)] != 1) == 1
            
            qk_design <- model.matrix(q_all[[ell]]$q_fits[[k]]@modelObj@model, data=df)  
            
            psi.betak <- diag(c(q_all[[ell]]$mod_regime_vhats[,k+1]*ones - q_all[[ell]]$unmod_regime_vhats[,k]*ones)) %*%
              qk_design
            psi.beta.ell <- cbind(psi.beta.ell, psi.betak) 

          } # end handling q_list length;
        } # end the if/else for existence of responders w/i feasible sets
        
      } # end the if/else of feasible sets, 
      
    } 
    psi.beta <- cbind(psi.beta, psi.beta.ell)
  }

  return(psi.beta)
}


#' Compute the negative derivative of the beta estimating equations with respect
#' to beta
#'
#' Computes \eqn{-\partial \Psi_{\beta} / \partial \beta}, the negative
#' derivative of the outcome regression estimating equations with respect to the
#' regression coefficients. This forms a block diagonal component of the An
#' matrix in the sandwich variance estimator.
#'
#' @param df A data frame containing the trial data, including treatment
#'   assignments, covariates, outcomes, and a \code{kappa} column.
#' @param q_all A list of fitted outcome regression objects for each regime,
#'   as returned by \code{\link{estimate_values}}.
#' @param regime_all A list of regime objects, each containing a \code{regime}
#'   matrix and a \code{regime_ind} indicator matrix.
#' @param feasible_sets_indicator A logical value indicating whether feasible sets
#'   are present in the trial design.
#' @param q_list A list of outcome regression model specifications (one per
#'   stage), used to determine the number of Q-function models at each stage
#'   (e.g., separate models for responders and non-responders).
#'
#' @return A block diagonal matrix (potentially sparse via \code{Matrix::bdiag})
#'   representing \eqn{-\partial \Psi_{\beta} / \partial \beta} across all
#'   regimes and stages.
#'
#' @noRd
ee_dpsi_beta <- function(df, q_all, regime_all, feasible_sets_indicator, q_list){
  for (ell in c(1:length(q_all))){
    K <- length(q_all[[ell]]$q_fits)
    
    for (k in c(1:K)){
      
      #### diagonals ####
      
      if (!feasible_sets_indicator){
        # multiple options for each stage every time
        ones <- (df['kappa']>=(k+1)) # those who have finished
        new <- (df['kappa']>=k) # those who can attain a prediction from model  
        
        ## diags
        qk_design <- model.matrix(q_all[[ell]]$q_fits[[k]]@modelObj@model, data=df[ones,])  # this is X
        betak_betak <- -t(qk_design) %*% qk_design
        
        ## handle matrix updates here
        
        
      } else{# begin when feasible sets occurs
        if (!(paste0("r", k) %in% colnames(df))){ # are we at a stage where
          # there are Rk columns? If no, (this part), then treat as normal
          # multiple options for each stage every time
          ones <- (df['kappa']>=(k+1)) # those who have finished
          new <- (df['kappa']>=k) # those who can attain a prediction from model  
          
          ## diags
          qk_design <- model.matrix(q_all[[ell]]$q_fits[[k]]@modelObj@model, data=df[ones,])  # this is X
          betak_betak <- -t(qk_design) %*% qk_design
          
        } else{ # there are rk at this stage
          # if this is true, then we need to determine if there are multiple 
          # models to fit. This occurs for interim analyses when we have some
          # individuals who have made it to stage k, had response determined,
          # but have not finished. 
          
          if (length(q_list[[k]]) > 1){
            ### start with those for rk==0 # non-responders ###
            ones <- (df['kappa']>=(k+1))*(df[paste0("r", k)] != 1) == 1
            new <- (df['kappa']>=k)*(df[paste0("r", k)] != 1) == 1
            
            ## diags
            qk_design_r0 <- model.matrix(q_all[[ell]]$q_fits[[k]][["r0"]]@modelObj@model, data=df[ones,])  # this is X
            betak_betak_r0 <- -t(qk_design_r0) %*% qk_design_r0
            
            ### now do responders rk == 1 ###
            # we fit on k+1 individuals, but only use the model hats on stage k
            ones <- (df['kappa']==(k+1))*(df[paste0("r", k)] == 1) == 1
            new <- (df['kappa']==k)*(df[paste0("r", k)] == 1) == 1
            
            ## diags
            qk_design_r1 <- model.matrix(q_all[[ell]]$q_fits[[k]][["r1"]]@modelObj@model, data=df[ones,])  # this is X
            betak_betak_r1 <- -t(qk_design_r1) %*% qk_design_r1
            
            # this should still be diagonal; hence off diagonals
            betak_betak <- Matrix::bdiag(betak_betak_r0, betak_betak_r1)
            
            
          } else{ # this section may occur if we are estimating AIPW, no interim
            ### start with those for rk==0 # non-responders ###
            ones <- (df['kappa']>=(k+1))*(df[paste0("r", k)] != 1) == 1
            new <- (df['kappa']>=k)*(df[paste0("r", k)] != 1) == 1
            
            ## diags
            qk_design <- model.matrix(q_all[[ell]]$q_fits[[k]]@modelObj@model, data=df[ones,])  # this is X
            betak_betak <- -t(qk_design) %*% qk_design
            
          } # end handling q_list length;
        } # end the if/else for existence of responders w/i feasible sets
        
      } # end the if/else of feasible sets, 
      
      
      #### off diagonals ####
      
      # in these cases, if we have feasible sets, we are interested in knowing 
      # if the response occurred a step above the current k and what was used 
      # for fitting the models 
      
      # this also should only occur if k != K
      if(k!=K){
        mod_df <- df
        # update the regime to match what is used to get response at the next stage
        mod_df[[paste0('a', k+1)]] <- regime_all[[ell]]$regime[,k+1]  
        
        if (!feasible_sets_indicator){
          # multiple options for each stage every time
          fit_ones <- (df['kappa']>=(k+1)) # those who have finished
          # then get the qkp1 design matrix, zeroed out for those who used Y to fit
          qkp1_design <- model.matrix(q_all[[ell]]$q_fits[[k+1]]@modelObj@model, data=mod_df[fit_ones,]) 
          qk_design <- model.matrix(q_all[[ell]]$q_fits[[k]]@modelObj@model, data=df[fit_ones,])
          
          betak_betakp1 <- t(qk_design) %*% qkp1_design
          
        } else{# begin when feasible sets occurs
          if (!(paste0("r", k+1) %in% colnames(df))){ # are we at a stage where
            # there are Rk columns? If no, (this part), then treat as normal
            # multiple options for each stage every time
            fit_ones <- (df['kappa']>=(k+1)) # those who have finished
            # then get the qkp1 design matrix, zeroed out for those who used Y to fit
            qkp1_design <- model.matrix(q_all[[ell]]$q_fits[[k+1]]@modelObj@model, data=mod_df[fit_ones,]) 
            qk_design <- model.matrix(q_all[[ell]]$q_fits[[k]]@modelObj@model, data=df[fit_ones,])
            
            betak_betakp1 <- t(qk_design) %*% qkp1_design
            
          } else{ # there are rk at this stage
            # if this is true, then we need to determine if there are multiple 
            # models to fit. This occurs for interim analyses when we have some
            # individuals who have made it to stage k, had response determined,
            # but have not finished. 
            
            if (length(q_list[[k+1]]) > 1){
              # this means that we have multiple functions to fit
              ### start with those for rk==0 # non-responders ###
              fit_ones <- (df['kappa']>=(k+1))*(df[paste0("r", k+1)] != 1) == 1
              qkp1_design_r0 <- model.matrix(q_all[[ell]]$q_fits[[k+1]][["r0"]]@modelObj@model, data=mod_df[fit_ones,]) 
              qk_design_r0 <- model.matrix(q_all[[ell]]$q_fits[[k]]@modelObj@model, data=df[fit_ones,])
              
              betak_betakp1_r0 <- t(qk_design_r0) %*% qkp1_design_r0
              
              
              ### now do responders rk == 1 ###
              # we fit on k+1 individuals, but only use the model hats on stage k
              fit_ones <- (df['kappa']==(k+1))*(df[paste0("r", k+1)] == 1) == 1
              qkp1_design_r1 <- model.matrix(q_all[[ell]]$q_fits[[k+1]][["r1"]]@modelObj@model, data=mod_df[fit_ones,]) 
              qk_design_r1 <- model.matrix(q_all[[ell]]$q_fits[[k]]@modelObj@model, data=df[fit_ones,])
              
              betak_betakp1_r1 <- t(qk_design_r1) %*% qkp1_design_r1
              
              
              betak_betakp1 <- cbind(betak_betakp1_r0, betak_betakp1_r1)
              # betak_betakp1 <- rbind(cbind(betak_betakp1_r0, matrix(0, nrow=nrow(betak_betakp1_r0), 
              #                                                 ncol=ncol(betak_betakp1_r1)) ),
              #                        cbind(matrix(0, nrow=nrow(betak_betakp1_r1), 
              #                                ncol=ncol(betak_betakp1_r0)), betak_betakp1_r1 ) )
              
            } else{ # this section may occur if we are estimating AIPW, no interim
              ### start with those for rk==0 # non-responders ###
              fit_ones <- (df['kappa']>=(k+1))*(df[paste0("r", k+1)] != 1) == 1
              # then get the qkp1 design matrix, zeroed out for those who used Y to fit
              qkp1_design <- model.matrix(q_all[[ell]]$q_fits[[k+1]]@modelObj@model, data=mod_df[fit_ones,]) 
              qk_design <- model.matrix(q_all[[ell]]$q_fits[[k]]@modelObj@model, data=df[fit_ones,])
              
              betak_betakp1 <- t(qk_design) %*% qkp1_design
              
            } # end handling q_list length;
          } # end the if/else for existence of responders w/i feasible sets
          
        } # end the if/else of feasible sets, 
        
      } # end if k!=K
      
      #### concat the matrices ####
      
      # we may have 2 matrices that occur, 
      
      if (k != K){
        dpsibetak <- cbind(betak_betak, betak_betakp1)
        trailingzeros <- dim(betak_betakp1)[2]
      } else{
        dpsibetak <- betak_betak
      }
      
      if (k==1){
        dpsi.beta.ell <- dpsibetak
      } else{
        # if k==K then we don't need the trailing zeros
        # we willl still need leadingzeros
        leadingzeros <- ncol(dpsi.beta.ell) - ncol(betak_betak)
        if (k==K){
          dpsi.beta.ell <- rbind(dpsi.beta.ell, 
                                 cbind(matrix(0, nrow=nrow(dpsibetak), 
                                              ncol=leadingzeros), 
                                       dpsibetak)) 
        } else{
          dpsi.beta.ell <- rbind(cbind(dpsi.beta.ell, matrix(0, nrow=nrow(dpsi.beta.ell), 
                                                             ncol=trailingzeros)), 
                                 cbind(matrix(0, nrow=nrow(dpsibetak), 
                                              ncol=leadingzeros), 
                                       dpsibetak)) 
        }
      }
      
    } # end for k in 1:K
    
    if (ell == 1) {
      dpsi.beta <- dpsi.beta.ell
    } else{
      dpsi.beta <- Matrix::bdiag(dpsi.beta, dpsi.beta.ell)
    }
    
  } # for loop for ells
  return(dpsi.beta)
}



