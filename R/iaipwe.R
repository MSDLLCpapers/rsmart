#' IAIPWE for K-stage SMARTs with up to 2 treatment options at each stage
#'
#' The IAIPWE was developed for arbitrary K stage SMART and subsumes both IPWE
#' and AIPWE. To implement the IPWE, q_list should be NULL and the t_s
#' should be set to the maximum available time of the outcome observed.
#' To implement the AIPWE, t_s should be set to the maximum available time of
#' the outcome observed.
#'
#' If times are not recorded for the study, "dummy" times can be used with the
#' analysis time set to be the maximum of the dummy times.
#'
#' @param df A data frame containing the data. It should include columns for
#'   the treatment assignments, response status, covariates, and outcomes.
#' @param pi_list A list of `modelObj` objects of length K specifying the
#'   propensity score model at each stage.
#' @param q_list A list of `modelObj` objects of length K specifying the
#'   outcome regression model at each stage, or `NULL` for IPW-only estimation.
#' @param regime_all A list of length equal to the number of regimes to be
#'   estimated. Each element is a list with two components: `regime`, an
#'   n x K matrix of treatment assignments each individual would receive under
#'   that regime, and `regime_ind`, an n x K indicator matrix of whether each
#'   individual was consistent with that regime at each stage.
#' @param feasible_sets_indicator A logical value indicating whether feasible
#'   sets are present in the trial design. If `TRUE`, some treatments are
#'   assigned deterministically based on response status; if `FALSE`, all
#'   stages have random treatment assignment.
#' @param t_s A numeric value indicating the time point at which the analysis
#'   should occur.
#' @param B A positive integer specifying the number of empirical bootstrap
#'   samples to use for variance estimation, or `NULL` (default) to use the
#'   asymptotic sandwich variance estimator.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{values}{A numeric vector of estimated regime values.}
#'     \item{se}{A numeric vector of standard errors for the regime values.}
#'     \item{covariance}{The estimated L x L covariance matrix of the regime
#'       value estimators.}
#'     \item{params}{A numeric vector of all estimated parameters.}
#'     \item{nus}{A list of estimated stage arrival probabilities, as returned
#'       by [get_nu()].}
#'     \item{q_all}{A list of fitted Q-function objects for each regime.}
#'     \item{regime_all}{The input regime list, returned for convenience.}
#'     \item{dfs}{A list of value term matrices for each regime.}
#'     \item{variance_choice}{A character string indicating the variance
#'       estimation method used (`"Asymptotic"` or `"Bootstrap"`).}
#'     \item{chi_square}{A list with `Statistic` (the chi-squared test
#'       statistic for equality of regime values) and `p_value`.}
#'   }
#'
#' @importFrom Matrix bdiag
#' @importFrom data.table rbindlist
#' @importFrom doFuture %dofuture%
#' @importFrom foreach foreach
#'
#' @export
#'

iaipwe <- function(df, pi_list, q_list, regime_all, feasible_sets_indicator, t_s,
                   B=NULL){
  # Check if the input data frame is empty
  if (nrow(df) == 0) {
    stop("Input data frame is empty.")
  }

  # Check if pi_list and q_list are of the same length
  if (!is.null(q_list) & length(pi_list) != length(q_list)) {
    stop("pi_list and q_list must have the same length.")
  }

  # Check if regimes is a valid vector
  if (!is.vector(regime_all) || length(regime_all) == 0) {
    stop("regime_all must be a non-empty vector.")
  }

  # Check if feasible_sets_indicator is a logical value
  if (!is.logical(feasible_sets_indicator)) {
    stop("feasible_sets_indicator must be a logical value.")
  }

  # Check if t_s is a numeric value
  if (!is.numeric(t_s) || length(t_s) != 1) {
    stop("t_s must be a single numeric value.")
  }

  # Check if B is NULL or a positive integer
  if (!(is.null(B) || (is.numeric(B) && length(B) == 1 && B > 1 && B %% 1 == 0))) {
    stop("B must be a positive integer greater than 1 or NULL.")
  }

  # Proceed with the function logic
  K <- length(pi_list)
  n <- nrow(df)

  #----------------------------------------------
  # estimate the values and q-functions as normal

  # estimate kappa
  kappa <- get_kappa(df, t_s, K)

  # if no one has finished do not run
  stopifnot(max(kappa) == K+1)

  # estimate nu
  df[,'kappa'] <- kappa
  nus <- get_nu(df, K)
  # estimate pis
  piFitted <- pi_fits(df, pi_list)
  pis <- piFitted[['ps']]
  p_fits <- piFitted[['p_fits']]
  # estimate values
  vHats<- estimate_values(df, q_list, regime_all, feasible_sets_indicator,
                          pis, nus)
  values <- vHats[['value']]
  dfs <- vHats[['df']]
  q_all <- vHats[['q_all']]

  # get all estimated parameters
  if (is.null(q_list)){
    estimatedParameters <- c(unlist(lapply(p_fits, function(x) x@fitObj$coefficients)),
                             unlist(nus$nu),
                             values
    )
  } else{
    estimatedParameters <- c(unlist(lapply(p_fits, function(x) x@fitObj$coefficients)),
                             unlist(nus$nu),
                             get_q_coefs(q_all),
                             values
    )
  }

  # ----------------------------------------------
  # If B is null, we will return the asymptotic variance
  if (is.null(B)){
    vartype <- "Asymptotic"
    # estimate the variance; these are already divided by n-p or n, respectively
    Bn <- get_bn(df, p_fits, nus, q_all, dfs, values, feasible_sets_indicator, q_list)
    An <- get_an(df, pis, p_fits, nus, q_all, values, regime_all, dfs, feasible_sets_indicator, q_list)

    # An is a sparse Matrix object, it may be singular.
    # In the event it is singular, take the ginv
    AnInv <- MASS::ginv(as.matrix(An))

    Vn <- AnInv %*% Bn %*% t(AnInv)
    # note that \widehat{\btheta} \sim N(0, Vn / n). Here n is actually nus$n
    # so the variance of \btheta is Vn / n
    L <- length(values)
    se <- sqrt(diag(Vn[(dim(Vn)[2] - L + 1) : dim(Vn)[2], (dim(Vn)[2] - L + 1) : dim(Vn)[2]]) / nus$ns)
  }

  # ----------------------------------------------
  # If B is not null, we will return the bootstrap variance
  else {
    vartype <- "Bootstrap"
    # foreach will return these items as a list; we use data.table to bind them
    boot_value_df <- foreach(
      i = seq_len(B),
      .errorhandling = "stop",
      .options.future = list(seed = TRUE)
    ) %dofuture% {
      # sample the rows of the data
      bootids <- sample(x = 1:nrow(df), size = nrow(df), replace = TRUE)
      # get boot dataset
      df_boot <- df[bootids,]
      # update regimes for the data set
      boot_regimes <- lapply(X = regime_all,
                             FUN = function(x) list("regime" = as.matrix(data.frame(x[[1]])[bootids,], rownames=FALSE),
                                                    "regime_ind" = as.matrix(data.frame(x[[2]])[bootids,], rownames=FALSE)))
      ## rename the rows
      rownames(df_boot) <- 1:nrow(df_boot)

      # estimate values for df_boot
      boot_nus <- get_nu(df_boot, K)
      # estimate pis
      boot_piFitted <- pi_fits(df_boot, pi_list)
      boot_pis <- boot_piFitted[['ps']]
      boot_p_fits <- boot_piFitted[['p_fits']]

      # estimate values
      boot_vHats <- estimate_values(df_boot, q_list, boot_regimes, feasible_sets_indicator,
                                    boot_pis, boot_nus)

      boot_values <- boot_vHats[['value']]
      as.data.frame(t(boot_values))
    }
    Vn <- cov(rbindlist(boot_value_df))
    se <- sqrt(diag(Vn))  # cov function already scales by B-1
  }

  L <- length(values)
  cont.mat <- cbind(diag(L-1), -1)
  endi <- dim(Vn)[2]
  starti <- endi - L + 1
  covvhat <- Vn[starti:endi, starti:endi]
  # only scale for the asymptotic error
  if (vartype == "Asymptotic") {
    covvhat <- covvhat / nus$ns
  }

  chisq <- t(cont.mat %*% values) %*%
    (solve(cont.mat %*% covvhat %*% t(cont.mat))) %*%
    (cont.mat %*% values)

  return(list('values' = values,
              'se' = se,
              'covariance' = covvhat,
              'params' = estimatedParameters,
              'nus' = nus,
              "q_all" = q_all,
              "regime_all" = regime_all,
              "dfs" = dfs,
              "variance_choice" = vartype,
              "chi_square" = list("Statistic" = chisq,
                                  "p_value" = 1-pchisq(q = chisq, df = as.integer(sum(svd(covvhat)$d > 1e-16))) ) ))
}
