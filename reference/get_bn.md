# Compute the Bn matrix for the sandwich variance estimator

Computes the Bn matrix, defined as \\n^{-1} \sum\_{i=1}^{n} \Psi_i
\Psi_i^T\\, the empirical variance of the estimating equations. Each row
of the individual-level estimating equation matrix corresponds to one
subject, with columns ordered as \\\pi_1, \ldots, \pi_K\\, \\\nu_1,
\ldots, \nu\_{K+1}\\, \\\beta^{(\ell)}\_1, \ldots, \beta^{(\ell)}\_K\\
for each regime \\\ell\\, and the value estimating equations \\V_1,
\ldots, V_L\\. Based on Section 7 of Boos and Stefanski for the robust
sandwich matrix.

## Usage

``` r
get_bn(df, p_fits, nus, q_all, dfs, values, feasible_sets_indicator, q_list)
```

## Arguments

- df:

  A data frame containing the trial data, including treatment
  assignments, covariates, outcomes, and a `kappa` column.

- p_fits:

  A list of fitted propensity score model objects (one per stage), each
  a `modelObj` fit object.

- nus:

  A list as returned by
  [`get_nu`](https://msdllcpapers.github.io/rsmart/reference/get_nu.md),
  containing the estimated stage probabilities `nu`, the sample size
  `ns`, and `nd`.

- q_all:

  A list of fitted outcome regression objects for each regime, as
  returned by
  [`estimate_values`](https://msdllcpapers.github.io/rsmart/reference/estimate_values.md).
  Can be an empty list when using IPW estimation.

- dfs:

  A list of data frames of value term components for each regime, as
  returned by
  [`estimate_values`](https://msdllcpapers.github.io/rsmart/reference/estimate_values.md).

- values:

  A numeric vector of estimated regime values.

- feasible_sets_indicator:

  A logical value indicating whether feasible sets are present in the
  trial design (i.e., some treatments are deterministic based on
  response status).

- q_list:

  A list of outcome regression model specifications (one per stage),
  used to determine the number of Q-function models at each stage (e.g.,
  separate models for responders and non-responders).

## Value

A numeric matrix representing the Bn component of the sandwich variance
estimator, with dimensions equal to the total number of estimated
parameters.
