#' Simulate treatment assignments for a single stage
#'
#' Generates treatment assignments for \code{n} individuals given a specified
#' number of possible treatments and a randomization probability function. This
#' is a utility function for constructing SMART simulations with flexible
#' randomization schemes. If an input data frame is provided, the treatment
#' assignments are appended as a new column \code{a<stage>}.
#'
#' @param n A positive integer. Number of individuals to assign treatments to.
#'   If \code{dat} is provided, \code{n} is inferred from \code{nrow(dat)} and
#'   this argument is ignored.
#' @param n_treatments A positive integer. Number of possible treatments.
#'   Treatments are coded as integers \code{0, 1, ..., n_treatments - 1}.
#' @param rand_prob_fn A function that takes \code{n} (the number of
#'   individuals), \code{n_treatments} (the number of treatment options), and
#'   \code{prob} (a numeric vector of probabilities) as arguments, and returns
#'   an integer vector of length \code{n} with treatment assignments coded as
#'   \code{0, 1, ..., n_treatments - 1}. Default is Bernoulli randomization
#'   using the probabilities specified by \code{prob}.
#' @param dat An optional data frame. If provided, the function appends the
#'   treatment assignments as a new column named \code{a<stage>} and returns
#'   the modified data frame. The data frame must not already contain a column
#'   with that name. Required when \code{stage > 1}.
#' @param stage A positive integer indicating the stage number. The treatment
#'   column will be named \code{paste0("a", stage)}. Default is 1. When
#'   \code{stage > 1}, \code{dat} must be provided and must contain columns
#'   \code{a1, ..., a<stage-1>}. Randomization is performed independently
#'   within each unique combination of prior treatments.
#' @param randomize_response One of \code{NULL} (default), \code{"Y"}, or
#'   \code{"N"}. If \code{"Y"}, the response indicator \code{r<stage>} is
#'   included in the grouping variable along with prior treatments
#'   \code{a1, ..., a<stage-1>} when randomizing. This allows different
#'   randomization within responder and non-responder subgroups.
#'   If \code{"N"}, responders (\code{r<stage> == 1}) are assigned treatment
#'   \code{0} deterministically, and only non-responders are randomized within
#'   prior treatment groups.
#'   Both \code{"Y"} and \code{"N"} require \code{stage > 1} and \code{dat} to
#'   contain a column named \code{r<stage>}. If \code{NULL}, response status is
#'   not used.
#' @param prob A numeric vector or a named list specifying randomization
#'   probabilities. If a numeric vector, it is used as the probability weights
#'   for all groups (must have length equal to \code{n_treatments}). If a named
#'   list, the names should correspond to the group levels (formed by
#'   \code{interaction()} of prior treatment columns and optionally the
#'   response column), and each element should be a numeric vector of length
#'   \code{n_treatments} giving the group-specific probabilities. Default is
#'   \code{NULL}, which uses equal probability across treatments.
#'
#' @return If \code{dat} is \code{NULL} (default), an integer vector of length
#'   \code{n} with treatment assignments coded as
#'   \code{0, 1, ..., n_treatments - 1}. If \code{dat} is provided, the input
#'   data frame with a new integer column \code{a<stage>} appended.
#'
#' @export
#'
#' @examples
#' # Default: Bernoulli randomization with 2 treatments (prob 0.5)
#' set.seed(1)
#' a <- sim_treatment(n = 100, n_treatments = 2)
#' table(a)
#'
#' # Unequal randomization: 70% to treatment 0, 30% to treatment 1
#' set.seed(1)
#' a <- sim_treatment(n = 100, n_treatments = 2, prob = c(0.7, 0.3))
#' table(a)
#'
#' # Three treatments with equal probability
#' set.seed(1)
#' a <- sim_treatment(n = 300, n_treatments = 3)
#' table(a)
#'
#' # With an input data frame
#' set.seed(1)
#' df <- data.frame(x1 = rnorm(100), x2 = rbinom(100, 1, 0.5))
#' df <- sim_treatment(n_treatments = 2, dat = df)
#' head(df)
#'
#' # With an input data frame at stage 2
#' set.seed(1)
#' df <- data.frame(x1 = rnorm(100), a1 = rbinom(100, 1, 0.5))
#' df <- sim_treatment(n_treatments = 2, dat = df, stage = 2)
#' head(df)
#'
#' \dontrun{
#' # Stage 2 with different randomization probabilities depending on a1:
#' # equal (50/50) if a1 == 0, unequal (70/30) if a1 == 1
#' # Using a named list for prob, where names match group levels
#' set.seed(1)
#' n <- 200
#' df <- data.frame(a1 = rbinom(n, 1, 0.5))
#' df <- sim_treatment(n_treatments = 2, dat = df, stage = 2,
#'                     prob = list("0" = c(0.5, 0.5), "1" = c(0.7, 0.3)))
#' table(df$a1, df$a2)
#' }
sim_treatment <- function(n = NULL, n_treatments,
                          rand_prob_fn = function(n, n_treatments, prob) {
                            sample(0:(n_treatments - 1), size = n,
                                   replace = TRUE, prob = prob)
                          },
                          dat = NULL,
                          stage = 1,
                          randomize_response = NULL,
                          prob = NULL) {
  # --- Input validation ---
  if (!is.numeric(stage) || stage < 1 || length(stage) != 1 || stage != floor(stage)) {
    stop("stage must be a positive integer.")
  }
  col_name <- paste0("a", stage)

  # Validate randomize_response

  if (!is.null(randomize_response)) {
    if (!(randomize_response %in% c("Y", "N"))) {
      stop("randomize_response must be NULL, 'Y', or 'N'.")
    }
    if (randomize_response %in% c("Y", "N") && stage == 1) {
      stop("randomize_response = 'Y' or 'N' requires stage > 1.")
    }
  }

  # If stage > 1, dat is required and must contain a1, ..., a(stage-1)
  if (stage > 1) {
    if (is.null(dat)) {
      stop("dat must be provided when stage > 1.")
    }
    prior_cols <- paste0("a", seq_len(stage - 1))
    missing_cols <- setdiff(prior_cols, colnames(dat))
    if (length(missing_cols) > 0) {
      stop(paste0("dat is missing prior treatment columns: ",
                  paste(missing_cols, collapse = ", "), "."))
    }
    # Check for response column if randomize_response == "Y" or "N"
    if (!is.null(randomize_response) && randomize_response %in% c("Y", "N")) {
      resp_col <- paste0("r", stage)
      if (!(resp_col %in% colnames(dat))) {
        stop(paste0("dat is missing response column '", resp_col,
                    "' required when randomize_response = '",
                    randomize_response, "'."))
      }
    }
  }

  if (!is.null(dat)) {
    if (!is.data.frame(dat)) {
      stop("dat must be a data frame.")
    }
    if (col_name %in% colnames(dat)) {
      stop(paste0("dat already contains a column named '", col_name, "'."))
    }
    n <- nrow(dat)
  }
  if (is.null(n) || !is.numeric(n) || n < 1 || length(n) != 1 || n != floor(n)) {
    stop("n must be a positive integer (or provide dat).")
  }
  if (!is.numeric(n_treatments) || n_treatments < 2 || length(n_treatments) != 1 ||
      n_treatments != floor(n_treatments)) {
    stop("n_treatments must be an integer >= 2.")
  }
  if (!is.function(rand_prob_fn)) {
    stop("rand_prob_fn must be a function.")
  }

  # Validate prob
  if (is.null(prob)) {
    prob <- rep(1 / n_treatments, n_treatments)
  }
  if (is.numeric(prob)) {
    if (length(prob) != n_treatments) {
      stop("prob vector must have length equal to n_treatments.")
    }
  } else if (is.list(prob)) {
    for (nm in names(prob)) {
      if (!is.numeric(prob[[nm]]) || length(prob[[nm]]) != n_treatments) {
        stop(paste0("Each element of prob list must be a numeric vector of ",
                    "length n_treatments. Problem with element '", nm, "'."))
      }
    }
  } else {
    stop("prob must be NULL, a numeric vector, or a named list.")
  }

  # Helper to resolve prob for a given group level
  resolve_prob <- function(grp) {
    if (is.list(prob)) {
      if (!(grp %in% names(prob))) {
        stop(paste0("prob list is missing an entry for group '", grp, "'."))
      }
      return(prob[[grp]])
    }
    return(prob)
  }

  # --- Generate treatment assignments ---
  if (stage > 1) {
    # Randomize within each unique combination of prior treatments a1, ..., a(stage-1)
    # and optionally response r<stage>
    prior_cols <- paste0("a", seq_len(stage - 1))
    group_cols <- prior_cols
    resp_col <- paste0("r", stage)

    if (!is.null(randomize_response) && randomize_response == "Y") {
      # Include response in grouping: randomize within (a1,...,a(stage-1), r<stage>)
      group_cols <- c(group_cols, resp_col)
    }

    treatments <- integer(n)

    if (!is.null(randomize_response) && randomize_response == "N") {
      # Responders receive treatment 0; only non-responders are randomized
      responders <- dat[, resp_col] == 1
      treatments[responders] <- 0L

      # Randomize non-responders within prior treatment groups
      non_resp_idx <- which(!responders)
      if (length(non_resp_idx) > 0) {
        groups_nr <- interaction(dat[non_resp_idx, group_cols, drop = FALSE], drop = TRUE)
        for (grp in levels(groups_nr)) {
          idx <- non_resp_idx[which(groups_nr == grp)]
          treatments[idx] <- rand_prob_fn(length(idx), n_treatments, resolve_prob(grp))
        }
      }
    } else {
      # Randomize all individuals within groups
      groups <- interaction(dat[, group_cols, drop = FALSE], drop = TRUE)
      for (grp in levels(groups)) {
        idx <- which(groups == grp)
        treatments[idx] <- rand_prob_fn(length(idx), n_treatments, resolve_prob(grp))
      }
    }
  } else {
    treatments <- rand_prob_fn(n, n_treatments, prob)
  }

  # --- Validate output ---
  if (length(treatments) != n) {
    stop("rand_prob_fn must return a vector of length n.")
  }
  if (!all(treatments %in% 0:(n_treatments - 1))) {
    stop("rand_prob_fn must return values in 0, 1, ..., n_treatments - 1.")
  }

  # --- Return ---
  if (!is.null(dat)) {
    dat[, col_name] <- as.integer(treatments)
    return(dat)
  }

  return(as.integer(treatments))
}
