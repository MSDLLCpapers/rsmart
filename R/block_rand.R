#' Create a blocked randomization function
#'
#' Returns a randomization function compatible with \code{\link{sim_treatment}}
#' that implements permuted block randomization. Within each block, every
#' treatment appears exactly \code{block_rep} times, and the order within blocks
#' is randomly permuted. This ensures approximate balance across treatments at
#' any point during enrollment.
#'
#' @param block_rep A positive integer. Number of times each treatment appears
#'   in every block. The block size is \code{block_rep * n_treatments}. Default
#'   is 2.
#'
#' @return A function with signature \code{function(n, n_treatments, prob)}
#'   suitable for use as \code{rand_prob_fn} in \code{\link{sim_treatment}}.
#'   The \code{prob} argument is accepted but ignored, since blocked
#'   randomization enforces equal allocation within each block.
#'
#' @export
#'
#' @examples
#' # Use blocked randomization with sim_treatment
#' set.seed(1)
#' a <- sim_treatment(n = 100, n_treatments = 2,
#'                    rand_prob_fn = block_rand(block_rep = 2))
#' table(a)
#'
#' # Three treatments with block size 6 (block_rep = 2)
#' set.seed(1)
#' a <- sim_treatment(n = 300, n_treatments = 3,
#'                    rand_prob_fn = block_rand(block_rep = 2))
#' table(a)
#'
#' # Stage 2 blocked randomization within prior treatment groups
#' set.seed(1)
#' df <- data.frame(a1 = c(rep(0, 50), rep(1, 50)))
#' df <- sim_treatment(n_treatments = 2, dat = df, stage = 2,
#'                     rand_prob_fn = block_rand(block_rep = 3))
#' table(df$a1, df$a2)
block_rand <- function(block_rep = 2) {
  if (!is.numeric(block_rep) || block_rep < 1 || length(block_rep) != 1 ||
      block_rep != floor(block_rep)) {
    stop("block_rep must be a positive integer.")
  }

  function(n, n_treatments, prob) {
    block_size <- block_rep * n_treatments
    n_blocks <- ceiling(n / block_size)
    allocation <- c(sapply(
      X = seq_len(n_blocks),
      FUN = function(x) {
        sample(rep(0:(n_treatments - 1), block_rep),
               size = block_size, replace = FALSE)
      }
    ))
    return(allocation[seq_len(n)])
  }
}
