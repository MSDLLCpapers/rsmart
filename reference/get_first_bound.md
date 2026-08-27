# Compute the first stopping boundary for a group sequential design

Determines the first analysis stopping boundary that controls the
familywise type I error rate at a specified level, adjusting for the
multiplicity of multiple treatment regimes. Uses an iterative search to
find the boundary value. This function works only for superiority
(one-sided) testing.

## Usage

``` r
get_first_bound(
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
  conducted. Values should be between 0 and 1.

- spend_fn:

  A character string specifying the alpha spending function. Either
  `"OF"` (O'Brien-Fleming) or `"Pocock"`.

- corr:

  A correlation matrix of Z-statistics at the first time point,
  accounting for multiple regimes. Default is a 1x1 identity matrix.

- test_type:

  A character string specifying the type of test to be performed. Either
  `"one-sided"` or `"two-sided"`. Default is `"one-sided"`.

- lambda:

  A numeric value for the initial step size used in the iterative
  boundary search. Default is 0.1.

- tol:

  A numeric value specifying the convergence tolerance for the type I
  error boundary. Default is 1e-6.

- max_iter:

  A positive integer specifying the maximum number of iterations for the
  boundary search. Default is 1000.

## Value

A list with the following components:

- bound:

  The computed stopping boundary for the first analysis.

- typeI:

  The achieved type I error rate at the boundary.

- convergence:

  A logical value indicating whether the algorithm converged.

- iters:

  The number of iterations used.
