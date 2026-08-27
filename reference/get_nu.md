# Estimate stage arrival probabilities (nu)

Estimates the probability that an individual has reached each stage
given they are enrolled in the trial. Requires the data frame to have a
`kappa` column indicating the stage reached. The returned list uses
indexing such that `nu[[k]]` corresponds to stage \\k\\, and `nu[[K+1]]`
is the probability that an individual has their final outcome observed.

## Usage

``` r
get_nu(df, K)
```

## Arguments

- df:

  A data frame containing a `kappa` column indicating the stage reached
  by each individual.

- K:

  An integer specifying the number of treatment stages (decision
  points).

## Value

A list with the following components:

- nu:

  A list of length \\K+1\\ where `nu[[k]]` is the estimated probability
  that an individual has reached stage \\k\\ given enrollment.

- ns:

  The total number of individuals enrolled in the trial (`kappa > 0`).

- nd:

  The proportion of individuals who have reached their last treatment
  stage and have their outcome observed.
