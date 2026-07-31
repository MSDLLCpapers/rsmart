#' Generate sample data for a two-stage SMART where responders are not 
#' re-randomized
#'
#' This function generates data for a two-stage Sequential Multiple Assignment
#' Randomized Trial (SMART) where responders are not re-randomized in the second
#' stage. The function allows for different value patterns and treatment assignments.
#'
#' @param n A positive integer. Number of individuals to generate data for.
#' @param s2 A positive numeric value. Variance of the error term. Default is
#'   100.
#' @param block_rep A positive integer. Number of times to duplicate treatments
#'   per block in the permuted block randomization. Default is 2.
#' @param r2p A numeric value between 0 and 1. Probability of being a
#'   responder in the second stage. Default is 0.5.
#'
#' @return A data frame with `n` rows and the following columns:
#'   \describe{
#'     \item{t1}{Study day of enrollment.}
#'     \item{t2}{Study day of second-stage randomization.}
#'     \item{t3}{Study day of outcome assessment.}
#'     \item{a1}{First-stage treatment assignment (0 or 1).}
#'     \item{r2}{Response indicator at stage 2 (0 = non-responder, 1 =
#'       responder).}
#'     \item{a2}{Second-stage treatment assignment (0 or 1 for non-responders,
#'       0 for responders).}
#'     \item{x11}{Baseline covariate (continuous, range 25--75).}
#'     \item{x12}{Baseline covariate (binary, 0 or 1).}
#'     \item{x21}{Stage 2 covariate (continuous, range 0--1).}
#'     \item{y}{Continuous outcome (higher is better).}
#'     \item{id}{Individual identifier (1 to `n`).}
#'   }
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- gen_no_trt_resp(n=400, s2=100, block_rep=2, r2p = 0.5)
#' head(dat)

gen_no_trt_resp <- function(n, s2=100, block_rep=2, r2p = 0.5){
  # Check if n is a non-negative integer
  if (!is.numeric(n) || n < 0 || length(n) != 1 || n != floor(n)) {
    stop("n must be a non-negative integer.")
  }
  # check that s2 is a positive numeric
  if (!is.numeric(s2) || s2 <= 0 || length(s2) != 1) {
    stop("s2 must be a positive numeric value.")
  }
  # check that block_rep is a non-negative integer
  if (!is.numeric(block_rep) || block_rep < 0 || length(block_rep) != 1 || block_rep != floor(block_rep)) {
    stop("block_rep must be a non-negative integer.")
  }
  # check that r2p is a numeric between 0 and 1
  if (!is.numeric(r2p) || r2p < 0 || r2p > 1 || length(r2p) != 1) {
    stop("r2p must be a numeric value between 0 and 1.")
  }

  # We generate the enrollment times of individuals and follow-up times as follows
  t1 <- round(runif(n=n, min=0, max=0.365*3), 3) * 1000  # enrollment over 2 years, in days
  t2 <- t1 + 0.1*1000 + round(runif(n=n, min=-1, max=1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  t3 <- t2 + 0.1*1000 + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10

  dat <- data.frame(t1, t2, t3)
  # order by t1 to match the allocation schedule when generated
  dat <- dat[order(dat$t1), ]

  # Stage 1: permuted block randomization with 2 treatments
  dat <- sim_treatment(n_treatments = 2, dat = dat, stage = 1,
                       rand_prob_fn = block_rand(block_rep = block_rep))

  # response status r2 is binomial with probability r2p
  dat[,"r2"] <- rbinom(n=n, size=1, prob=r2p)

  # order by t2 to match the allocation schedule for stage 2
  dat <- dat[order(dat[,"t2"]),]

  # Stage 2: non-responders are re-randomized within a1 groups,
  # responders receive a2 = 0
  dat <- sim_treatment(n_treatments = 2, dat = dat, stage = 2,
                       rand_prob_fn = block_rand(block_rep = block_rep),
                       randomize_response = "N")

  # now we generate some covariates for the data to make the IAIPWE useful
  dat[,"x11"] <- round(runif(n=n, min=25, max=75)) # mean 50 # age
  dat[,"x12"] <- rbinom(n=n, size=1, prob=0.5) # mean 0.5 # gender
  dat[,"x21"] <- round(runif(n=n, min=0, max=1),2) # mean 0.5 # adherence

  # now we generate the errors
  ei <- rnorm(n=n, mean=0, sd=sqrt(s2))

  # we create the outcomes as a linear combination of the treatments, covariates, and errors
  # here, receiving a1=0 is better, responders also have a higher outcome, and second stage a2=0 is improved
  # baseline improvement for all regimes is 10
  # the covariates contribute an average of 30 to the outcome on top of the 10
  # the treatments contribute an additional 2 for participants who get a1=1 and 5 if a2 = 1
  # and an additional 3 if both a1 and a2 are 1
  # responders have an increase of 8, which matches that of the a1 and a2 = 1 group

  dat[,"y"] <- 10 + dat[,"a1"] * 2 + dat[,"a2"] * 5  +
    (dat[,"a1"] * dat[,"a2"]) * 3 +
    dat[,"r2"] * 8 +
    dat[,"x11"] * 0.5 + dat[,"x12"] * 5 + dat[,"x21"] * 5 +
    ei

  # the expected outcomes are as follows:
  # A1 | R2 | A2 | E(Y)
  #  0 |  0 |  0 |  10 + 0 + 0 + 0 + 0 + 25 + 2.5 + 2.5 = 40
  #  0 |  0 |  1 |  10 + 0 + 5 + 0 + 0 + 25 + 2.5 + 2.5 = 45
  #  0 |  1 |  0 |  10 + 0 + 0 + 0 + 8 + 25 + 2.5 + 2.5 = 48
  #  -------------------------------------------------------
  #  1 |  0 |  0 |  10 + 2 + 0 + 0 + 0 + 25 + 2.5 + 2.5 = 42
  #  1 |  0 |  1 |  10 + 2 + 5 + 3 + 0 + 25 + 2.5 + 2.5 = 50
  #  1 |  1 |  0 |  10 + 2 + 0 + 0 + 8 + 25 + 2.5 + 2.5 = 48
  #  -------------------------------------------------------
  # So Regimes have values 44, 46.5, 45, and 49.

  dat <- dat[order(dat$t1), ]

  # add an id for each individual, it should be numbered 1:n
  dat[,"id"] <- 1:n

  return(dat)
}
