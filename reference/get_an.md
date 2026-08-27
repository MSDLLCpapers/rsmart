# Compute the An matrix for the sandwich variance estimator

Computes the An matrix, defined as \\-\partial \Psi / \partial \theta\\,
where \\\Psi\\ is the stacked estimating equations vector. The
estimating equations are ordered to match Bn: \\\pi_1, \ldots, \pi_K\\,
\\\nu_1, \ldots, \nu\_{K+1}\\, \\\beta^{(\ell)}\_1, \ldots,
\beta^{(\ell)}\_K\\ for each regime \\\ell\\, and the value estimating
equations \\V_1, \ldots, V_L\\. Based on Section 7 of Boos and Stefanski
for the robust sandwich matrix.

## Usage

``` r
get_an(
  df,
  pis,
  p_fits,
  nus,
  q_all,
  values,
  regime_all,
  dfs,
  feasible_sets_indicator,
  q_list
)
```

## Arguments

- df:

  A data frame containing the trial data, including treatment
  assignments, covariates, outcomes, and a `kappa` column.

- pis:

  A data frame of estimated propensity scores with columns
  `pi1, ..., piK`, one column per stage.

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

- values:

  A numeric vector of estimated regime values.

- regime_all:

  A list of regime objects, each containing a `regime` matrix and a
  `regime_ind` indicator matrix.

- dfs:

  A list of data frames of value term components for each regime, as
  returned by
  [`estimate_values`](https://msdllcpapers.github.io/rsmart/reference/estimate_values.md).

- feasible_sets_indicator:

  A logical value indicating whether feasible sets are present in the
  trial design (i.e., some treatments are deterministic based on
  response status).

- q_list:

  A list of outcome regression model specifications (one per stage),
  used to determine the number of Q-function models at each stage (e.g.,
  separate models for responders and non-responders).

## Value

A numeric matrix representing the An component of the sandwich variance
estimator, with dimensions equal to the total number of estimated
parameters.
