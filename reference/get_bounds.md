# Compute all stopping boundaries for a two-analysis group sequential design

Wrapper function that computes the stopping boundaries for both the
first and second analyses of a group sequential design with multiple
treatment regimes. Calls
[`get_first_bound`](https://msdllcpapers.github.io/rsmart/reference/get_first_bound.md)
and
[`get_next_bound`](https://msdllcpapers.github.io/rsmart/reference/get_next_bound.md)
sequentially.

## Usage

``` r
get_bounds(
  alpha = 0.05,
  inf_frac = c(0.5, 1),
  spend_fn = "OF",
  corr = diag(x = 1, nrow = 1, ncol = 1),
  test_type = "one-sided",
  lambda = 0.1,
  tol = 1e-06,
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

  A correlation matrix of Z-statistics at analysis time s. Should have
  dimension \\L \times L\\.

- test_type:

  A character string specifying the type of test to be performed. Either
  `"one-sided"` or `"two-sided"`. Default is `"one-sided"`.

- lambda:

  A numeric value for the initial step size used in the iterative
  boundary search. Default is 0.1.

- tol:

  A numeric value specifying the convergence tolerance. Default is 1e-6.

- max_iter:

  A positive integer specifying the maximum number of iterations.
  Default is 1000.

## Value

A list with the following components:

- bounds / bound:

  A numeric value (single analysis) or numeric vector (sequential) of
  stopping boundaries. The key is `"bounds"` for a single analysis and
  `"bound"` for multiple analyses.

- spending:

  A numeric value or vector of cumulative alpha spent at each analysis.

- convergence:

  A logical value or vector indicating whether the algorithm converged
  at each analysis.

- iterations:

  An integer value or vector of the number of iterations used at each
  analysis.

- alpha:

  The input type I error rate.

- inf_frac:

  The input information fractions.

- spend_fn:

  The input spending function name.

- corr:

  The input correlation matrix.

- test_type:

  The input test type.
