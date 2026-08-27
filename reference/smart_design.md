# Compute stopping boundaries and sample size for a group sequential SMART

Calculates the stopping boundaries and required sample size for a group
sequential design with multiple treatment regimes embedded in a SMART.
The boundaries are found first to control the type I error rate, then
the sample size is determined to achieve the desired power. This
function wraps
[`get_bounds`](https://msdllcpapers.github.io/rsmart/reference/get_bounds.md)
and
[`get_sample_size`](https://msdllcpapers.github.io/rsmart/reference/get_sample_size.md).

## Usage

``` r
smart_design(
  test_type = "one-sided",
  comp_type = "fixed.control",
  alpha = 0.05,
  beta = 0.1,
  delta = NULL,
  n_init = 400,
  inf_frac = 1,
  spend_fn = "OF",
  corr = NULL,
  variances = NULL,
  lambdaB = 0.1,
  lambdaC = 20,
  tol = 1e-06,
  max_iterB = 1000,
  max_iterC = 100,
  seed = 1,
  Bb = 100001,
  Bc = 10001,
  mu = NULL
)
```

## Arguments

- test_type:

  A character string specifying the type of hypothesis test. A test type
  of either one-sided or two-sided. Two-sided testing must be symmetric.
  This input is ignored for a global Chi-squared test.

- comp_type:

  A character string specifying what comparison will be tested.
  Character string for either "fixed.control" or "global.chi.sq".

- alpha:

  A numeric value between 0 and 1 specifying the overall familywise type
  I error rate. Default is `0.05`.

- beta:

  A numeric value between 0 and 1 specifying the type II error rate.
  Power is `1 - beta`. Default is `0.1` for power of 0.9.

- delta:

  A numeric vector of length equal to the number of treatment regimes
  specifying the alternative differences between each regime's value and
  the null (or control) value.

- n_init:

  A positive integer specifying the initial total sample size to begin
  the iterative sample size search.

- inf_frac:

  A numeric vector of information fractions indicating when analyses are
  conducted. Values should be between 0 and 1, with the last element
  equal to 1. Default is `1` (single analysis). For example, `c(0.5, 1)`
  specifies an interim analysis at 50\\ information and a final
  analysis.

- spend_fn:

  A character string specifying the alpha spending function. Either
  `"OF"` (O'Brien-Fleming, default) or `"Pocock"`.

- corr:

  A correlation matrix of Z-statistics across all regimes and analyses.
  The dimension should be \\L \times k\\ by \\L \times k\\, where \\L\\
  is the number of treatment regimes. This matrix accounts for both
  within-analysis correlations between regimes and across-analysis
  correlations due to shared participants.

- variances:

  A numeric vector of length equal to the number of treatment regimes
  giving the variance of each value estimator scaled by sample size,
  i.e., \\N \times \mathrm{Var}(\hat{V}\_d)\\.

- lambdaB:

  A numeric value for the initial step size used in the iterative
  boundary search. Default is 0.1.

- lambdaC:

  A numeric value for the initial step size used in the iterative sample
  size search. Default is 20.

- tol:

  A numeric value specifying the convergence tolerance. Default is 1e-6.

- max_iterB:

  A positive integer specifying the maximum number of iterations for the
  boundary search. Default is 1000.

- max_iterC:

  A positive integer specifying the maximum number of iterations for the
  sample size search. Default is 100.

- seed:

  An integer seed for reproducibility of the Monte Carlo simulation.
  Default is 1.

- Bb:

  A positive integer specifying the number of Monte Carlo samples for
  the boundary computation. Default is 100001.

- Bc:

  A positive integer specifying the number of Monte Carlo samples for
  the sample size computation. Default is 10001.

- mu:

  A parameter vector for the multivariate normal distribution of the
  treatment effect mean, used only for the chi-square distribution in
  the case of a non-centrality assumption, though unlikely to be needed.

## Value

A list with the following components:

- boundaries:

  The output from
  [`get_bounds`](https://msdllcpapers.github.io/rsmart/reference/get_bounds.md)
  (when `comp_type = "fixed.control"`) or
  [`get_bounds_chi`](https://msdllcpapers.github.io/rsmart/reference/get_bounds_chi.md)
  (when `comp_type = "global.chi.sq"`), containing the stopping
  boundaries and associated metadata. See those functions for details.

- sample.size:

  The output from
  [`get_sample_size`](https://msdllcpapers.github.io/rsmart/reference/get_sample_size.md)
  or
  [`get_sample_size_chi`](https://msdllcpapers.github.io/rsmart/reference/get_sample_size_chi.md),
  containing the required sample sizes, achieved power, cumulative
  rejection probabilities, and echoed input parameters. If `delta` is
  `NULL` (no alternative specified), the sample size is not computed and
  this element is `NULL`.

- inputs:

  A list echoing all input parameters passed to `smart_design`.
