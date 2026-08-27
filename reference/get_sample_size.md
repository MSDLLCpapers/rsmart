# Determine sample size for a group sequential SMART design

For specified operating characteristics, iteratively increases the
sample size until the desired power is achieved. Assumes the information
fraction remains unchanged as the sample size increases, which is
reasonable for small changes or at the design stage.

## Usage

``` r
get_sample_size(
  variances,
  beta,
  delta,
  bounds,
  n_init = 100,
  corr = bounds$corr,
  inf_frac = bounds$inf_frac,
  n_split = NULL
)
```

## Arguments

- variances:

  A numeric vector of length \\L\\ of variances of the value estimators,
  i.e., \\\sqrt{N} \times \mathrm{Cov}(\hat{\theta})\\. These should
  reflect the population variances rather than sample variance or
  standard errors.

- beta:

  A numeric value between 0 and 1 specifying the type II error rate.
  Power is `1 - beta`.

- delta:

  A numeric vector of length \\L\\ of differences between regime values
  and the null value (or control arm).

- bounds:

  A numeric vector of length \\S\\ with the stopping boundaries for
  analyses \\1, \ldots, S\\, or the boundaries from
  [`get_bounds`](https://msdllcpapers.github.io/rsmart/reference/get_bounds.md)
  function

- n_init:

  A positive integer specifying the initial total trial sample size
  \\N\\ to begin the search.

- corr:

  A correlation matrix of dimension \\L \times L\\ between regime value
  estimators across all analyses.

- inf_frac:

  A numeric vector of information fractions indicating when analyses are
  conducted.

- n_split:

  A numeric vector indicating the proportion of the sample size at the
  analysis times s=1,...,S. If no argument is given, assumes the split
  is proportional to the information available.

## Value

A list with the following components:

- N:

  A numeric vector of sample sizes at each analysis.

- power:

  The achieved power at the final sample size.

- prop_rej:

  A numeric vector of cumulative rejection probabilities at each
  analysis.

- variances:

  The input variances.

- beta:

  The input type II error rate.

- delta:

  The input alternative differences.

- bounds:

  The input stopping boundaries.

- n_init:

  The input initial sample size.

- corr:

  The input correlation matrix.

- inf_frac:

  The input information fractions.

- n_split:

  The input sample size split proportions.
