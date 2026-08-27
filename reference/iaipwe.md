# IAIPWE for K-stage SMARTs with up to 2 treatment options at each stage

The IAIPWE was developed for arbitrary K stage SMART and subsumes both
IPWE and AIPWE. To implement the IPWE, q_list should be NULL and the t_s
should be set to the maximum available time of the outcome observed. To
implement the AIPWE, t_s should be set to the maximum available time of
the outcome observed.

## Usage

``` r
iaipwe(df, pi_list, q_list, regime_all, feasible_sets_indicator, t_s, B = NULL)
```

## Arguments

- df:

  A data frame containing the data. It should include columns for the
  treatment assignments, response status, covariates, and outcomes.

- pi_list:

  A list of `modelObj` objects of length K specifying the propensity
  score model at each stage.

- q_list:

  A list of `modelObj` objects of length K specifying the outcome
  regression model at each stage, or `NULL` for IPW-only estimation.

- regime_all:

  A list of length equal to the number of regimes to be estimated. Each
  element is a list with two components: `regime`, an n x K matrix of
  treatment assignments each individual would receive under that regime,
  and `regime_ind`, an n x K indicator matrix of whether each individual
  was consistent with that regime at each stage.

- feasible_sets_indicator:

  A logical value indicating whether feasible sets are present in the
  trial design. If `TRUE`, some treatments are assigned
  deterministically based on response status; if `FALSE`, all stages
  have random treatment assignment.

- t_s:

  A numeric value indicating the time point at which the analysis should
  occur.

- B:

  A positive integer specifying the number of empirical bootstrap
  samples to use for variance estimation, or `NULL` (default) to use the
  asymptotic sandwich variance estimator.

## Value

A list with the following components:

- values:

  A numeric vector of estimated regime values.

- se:

  A numeric vector of standard errors for the regime values.

- covariance:

  The estimated L x L covariance matrix of the regime value estimators.

- params:

  A numeric vector of all estimated parameters.

- nus:

  A list of estimated stage arrival probabilities, as returned by
  [`get_nu()`](https://msdllcpapers.github.io/rsmart/reference/get_nu.md).

- q_all:

  A list of fitted Q-function objects for each regime.

- regime_all:

  The input regime list, returned for convenience.

- dfs:

  A list of value term matrices for each regime.

- variance_choice:

  A character string indicating the variance estimation method used
  (`"Asymptotic"` or `"Bootstrap"`).

- chi_square:

  A list with `Statistic` (the chi-squared test statistic for equality
  of regime values) and `p_value`.

## Details

If times are not recorded for the study, "dummy" times can be used with
the analysis time set to be the maximum of the dummy times.
