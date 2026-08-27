# Fit outcome regression (Q-function) models across all stages for a single regime

Fits Q-function models backwards from the last stage to the first for a
single treatment regime. At each stage, both regime-modified and
unmodified predicted values are computed. Handles feasible sets (where
responders may not be re-randomized) and interim analyses.

## Usage

``` r
get_q_fits(df, q_list, regime, feasible_sets_indicator = FALSE)
```

## Arguments

- df:

  A data frame containing the trial data, including treatment
  assignments, covariates, outcomes (`y`), and a `kappa` column.

- q_list:

  A list of outcome regression model specifications (one per stage).
  Each element is either a `modelObj` object or a named list with
  elements `"r0"` and `"r1"` for non-responders and responders,
  respectively.

- regime:

  A matrix of treatment assignments under the regime being evaluated,
  with columns corresponding to stages.

- feasible_sets_indicator:

  A logical value indicating whether feasible sets are present. Default
  is `FALSE`. If individuals are not re-randomized (i.e. for some stage
  a response status prevents individuals from receiving a random
  treatment), then feasible_sets_indicator should be set to TRUE.

## Value

A list with the following components:

- q_fits:

  A list of fitted `modelObj` objects, one per stage. When feasible sets
  with multiple models are used, the element is a named list with `"r0"`
  and `"r1"`.

- mod_regime_vhats:

  A matrix of regime-modified predicted values with columns
  `q1, ..., q_{K+1}`.

- unmod_regime_vhats:

  A matrix of unmodified predicted values with columns
  `q1_nochange, ..., q_{K+1}_nochange`.
