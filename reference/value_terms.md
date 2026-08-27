# Compute the individual-level value terms for a single regime

Computes the \\2K+1\\ coarsening-level value terms for each individual
under a single treatment regime. These include augmentation terms for
levels \\r = 1, \ldots, 2K\\ and the IPW term for \\R = \infty\\ (i.e.,
\\R = 2K+1\\). When outcome regression models (`qs`) are provided, the
augmented terms are computed; otherwise only the IPW term is non-zero.

## Usage

``` r
value_terms(df, regime_ind, pis, qs, nus)
```

## Arguments

- df:

  A data frame containing the trial data, including treatment
  assignments, covariates, outcomes (`y`), and a `kappa` column.

- regime_ind:

  A matrix of regime consistency indicators for a single regime, with
  columns indicating whether each individual followed the regime at each
  stage.

- pis:

  A data frame of estimated propensity scores with columns
  `pi1, ..., piK`.

- qs:

  A list containing outcome regression fitted values for a single regime
  (as returned by
  [`get_q_fits`](https://msdllcpapers.github.io/rsmart/reference/get_q_fits.md)),
  or `NULL` for IPW estimation.

- nus:

  A list as returned by
  [`get_nu`](https://msdllcpapers.github.io/rsmart/reference/get_nu.md),
  containing the estimated stage probabilities.

## Value

A numeric matrix with `nrow(df)` rows and \\2K+1\\ columns, where each
column corresponds to a coarsening-level term in the value estimator.
