# Fit propensity score models for all stages

Fits propensity score models at each stage of the SMART, using only
individuals who have reached that stage. Returns estimated propensity
scores and fitted model objects for all stages.

## Usage

``` r
pi_fits(df, p_list)
```

## Arguments

- df:

  A data frame containing the trial data, including treatment
  assignments (`a1, a2, ...`), covariates, and a `kappa` column.

- p_list:

  A list of `modelObj` objects specifying the propensity score model at
  each stage.

## Value

A list with the following components:

- ps:

  A data frame of estimated propensity scores with columns
  `pi1, ..., piK`. Individuals who have not reached a stage are assigned
  a value of 99.

- p_fits:

  A list of fitted `modelObj` objects, one per stage.
