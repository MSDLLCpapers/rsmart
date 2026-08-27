# Compute the kappa (stage reached) for each individual

Determines the total number of stages each individual has reached by a
specified analysis time `t_s`, based on their arrival times
`t1, ..., t_{K+1}`.

## Usage

``` r
get_kappa(df, t_s, K)
```

## Arguments

- df:

  A data frame containing columns `t1, t2, ..., t_{K+1}` representing
  the times at which each individual reaches each stage.

- t_s:

  A numeric value specifying the analysis time point.

- K:

  An integer specifying the number of treatment stages (decision
  points).

## Value

An integer vector of length `nrow(df)` indicating the number of stages
each individual has reached by time `t_s`. Values range from 0 (not yet
enrolled) to \\K+1\\ (outcome observed).
