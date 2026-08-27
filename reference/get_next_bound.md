# Compute subsequent stopping boundaries for a group sequential design

Given the boundary from the first analysis, determines the stopping
boundary for the \\s\\-th analysis that controls the overall type I
error rate. Uses an iterative search with the joint distribution of test
statistics across analyses.

## Usage

``` r
get_next_bound(
  alpha = 0.05,
  inf_frac = c(0.5, 1),
  spend_fn = "OF",
  corr = diag(x = 1, nrow = 1, ncol = 1),
  test_type = "one-sided",
  lambda = 0.1,
  tol = 1e-06,
  prev_bound = NULL,
  s = 2,
  max_iter = 1000
)
```

## Arguments

- alpha:

  A numeric value specifying the overall type I error rate to control.
  Default is 0.05.

- inf_frac:

  A numeric vector of information fractions indicating when analyses are
  conducted.

- spend_fn:

  A character string specifying the alpha spending function. Either
  `"OF"` (O'Brien-Fleming) or `"Pocock"`.

- corr:

  A correlation matrix of Z-statistics across all analyses and regimes.
  Should have dimension \\SL \times SL\\ by \\SL \times SL\\.

- test_type:

  A character string specifying the type of test to be performed. Either
  `"one-sided"` or `"two-sided"`. Default is `"one-sided"`.

- lambda:

  A numeric value for the initial step size used in the iterative
  boundary search. Default is 0.1.

- tol:

  A numeric value specifying the convergence tolerance. Default is 1e-6.

- prev_bound:

  A numeric vector of previously computed boundaries from analyses \\1,
  \ldots, s-1\\, with dimension \\(s-1) \times L\\.

- s:

  An integer indicating the current analysis number. Default is 2.

- max_iter:

  A positive integer specifying the maximum number of iterations.
  Default is 1000.

## Value

A list with the following components:

- bound:

  The computed stopping boundary for analysis `s`.

- typeI:

  The achieved cumulative type I error rate through analysis `s`.

- convergence:

  A logical value indicating whether the algorithm converged.

- iters:

  The number of iterations used.
