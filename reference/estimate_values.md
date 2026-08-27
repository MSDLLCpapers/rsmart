# Estimate values for all treatment regimes

Wraps the outcome regression and value term functions to estimate the
value of all regimes in the provided list. For each regime, Q-functions
are fitted (if specified), augmentation and IPW terms are computed, and
the regime value is estimated.

## Usage

``` r
estimate_values(df, q_list, regime_all, feasible_sets_indicator, pis, nus)
```

## Arguments

- df:

  A data frame containing the trial data, including treatment
  assignments, covariates, outcomes, and a `kappa` column.

- q_list:

  A list of outcome regression model specifications (one per stage), or
  `NULL` for IPW estimation.

- regime_all:

  A list of regime objects, each containing a `regime` matrix and a
  `regime_ind` indicator matrix.

- feasible_sets_indicator:

  A logical value indicating whether feasible sets are present in the
  trial design.

- pis:

  A data frame of estimated propensity scores with columns
  `pi1, ..., piK`.

- nus:

  A list as returned by
  [`get_nu`](https://msdllcpapers.github.io/rsmart/reference/get_nu.md),
  containing the estimated stage probabilities.

## Value

A list with the following components:

- value:

  A numeric vector of estimated regime values.

- df:

  A list of matrices of value term components, augmentation terms 1
  through 2K then the IPW term for each regime.

- q_all:

  A list of fitted Q-function objects for each regime, or an empty list
  if `q_list` is `NULL`.
