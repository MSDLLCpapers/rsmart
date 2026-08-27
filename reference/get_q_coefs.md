# Extract coefficients from all fitted Q-function models

Extracts and concatenates the regression coefficients from all fitted
Q-function models across all regimes and stages.

## Usage

``` r
get_q_coefs(q_all)
```

## Arguments

- q_all:

  A list of fitted outcome regression objects for each regime, as
  returned by
  [`estimate_values`](https://msdllcpapers.github.io/rsmart/reference/estimate_values.md).

## Value

A numeric vector of all Q-function regression coefficients, ordered by
regime and then by stage.
