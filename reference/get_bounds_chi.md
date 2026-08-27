# Compute the first chi-squared stopping boundary for a group sequential design

Determines the first analysis stopping boundary based on a chi-squared
global test statistic. This is used when the comparison type is a global
chi-squared test (testing whether any regime differs from the others)
rather than individual regime comparisons against a fixed control.

## Usage

``` r
get_bounds_chi(
  alpha = 0.05,
  inf_frac = c(0.5, 1),
  spend_fn = "OF",
  corr = diag(x = 1, nrow = 1, ncol = 1),
  mu = NULL,
  B = 1000001,
  seed = 1
)
```

## Arguments

- alpha:

  A numeric value specifying the overall type I error rate to control.
  Default is 0.05.

- inf_frac:

  A numeric vector of information fractions indicating when analyses are
  conducted. Values should be between 0 and 1. The length determines the
  number of planned analyses \\S\\. Default is `c(0.5, 1)`.

- spend_fn:

  A character string specifying the alpha spending function. Currently
  used to adjust the `iota` scaling factors for each analysis boundary.
  For `"OF"` (O'Brien-Fleming), set `iota` to the information fractions.
  For `"Pocock"`, use `iota = rep(1, S)`. Default is `"OF"`.

- corr:

  A correlation matrix of the expected correlation between the
  Z-statistics of regimes against a fixed value. Should have dimension
  \\L \times L\\ by \\L \times L\\, where \\L\\ is the number of
  regimes. For `inf_frac` with length greater than one, the correlation
  structure across multiple analyses is computed.

- mu:

  An optional numeric vector of means for the multivariate normal
  distribution for all regimes at a single analysis used in the Monte
  Carlo simulation. Default is `NULL`, which uses a zero vector (null
  hypothesis).

- B:

  A positive integer specifying the number of Monte Carlo samples.
  Default is 1000001.

- seed:

  An integer seed for reproducibility of the Monte Carlo simulation.
  Default is 1.

## Value

A list with the following components:

- bound:

  A numeric vector of chi-squared stopping boundaries, one per analysis.

- spending:

  A numeric vector of observed cumulative alpha spent at each analysis.

- typeI:

  The overall achieved type I error rate across all analyses.

- convergence:

  A logical value; always `TRUE` for Monte Carlo based computation.

- iters:

  The number of Monte Carlo samples used (`B`).

- dfchi:

  The degrees of freedom of the chi-squared statistic, equal to the rank
  of \\C \Sigma C'\\.

- alpha:

  The input type I error rate.

- inf_frac:

  The input information fractions.

- spend_fn:

  The input spending function name.

- corr:

  The input correlation matrix.

- mu:

  The input mean vector.

- B:

  The input number of Monte Carlo samples.

- seed:

  The input random seed.

## Details

The function uses Monte Carlo simulation to estimate the boundary that
controls the familywise type I error rate at a specified level. A
contrast matrix \\C\\ is constructed to form \\L-1\\ linearly
independent comparisons among \\L\\ regimes, and the chi-squared
statistic is computed as \\(CZ)' (C \Sigma C')^{-1} (CZ)\\.
