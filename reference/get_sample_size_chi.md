# Determine sample size for a chi-squared global test in a group sequential SMART design

For specified operating characteristics, iteratively increases the
sample size until the desired power is achieved using a chi-squared
global test statistic. The chi-squared test assesses whether any
treatment regime differs from the others, rather than comparing
individual regimes against a fixed control.

## Usage

``` r
get_sample_size_chi(
  variances,
  beta,
  delta,
  bounds,
  n_init = 100,
  corr = bounds$corr,
  inf_frac = bounds$inf_frac,
  n_split = NULL,
  seed = 2,
  B = 10001,
  lambda = 20,
  max_iter = 100
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

  A numeric vector of length \\L\\ of change in the regime values under
  the alternative hypothesis.

- bounds:

  A numeric vector of chi-squared stopping boundaries for analyses \\1,
  \ldots, S\\, or the output list from
  [`get_bounds_chi`](https://msdllcpapers.github.io/rsmart/reference/get_bounds_chi.md).

- n_init:

  A positive integer specifying the initial total trial sample size
  \\N\\ to begin the search. Default is 100.

- corr:

  A correlation matrix of dimension \\L \times L\\ between regime value
  estimators. Defaults to the correlation from the `bounds` list if
  provided.

- inf_frac:

  A numeric vector of information fractions indicating when analyses are
  conducted. Defaults to the information fractions from the `bounds`
  list if provided.

- n_split:

  A numeric vector indicating the proportion of the sample size at the
  analysis times \\s = 1, \ldots, S\\. If `NULL` (default), assumes the
  split is proportional to the information fractions.

- seed:

  An integer seed for reproducibility of the Monte Carlo simulation.
  Default is 2.

- B:

  A positive integer specifying the number of Monte Carlo samples for
  power estimation. Default is 10001.

- lambda:

  A positive numeric value specifying the incremental sample size
  increases when raising the sample size to find the required power via
  simulation. After this is found, iteration from \\n-lambda\\ to
  \\n+lambda\\ will be done to find the correct exact sample size.

- max_iter:

  The maximum sample sizes n to evaluate.

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

- seed:

  The input random seed.

- B:

  The input number of Monte Carlo samples.

- lambda:

  The input step size for the sample size search.

## Details

The function uses Monte Carlo simulation to estimate power at each
candidate sample size. A contrast matrix \\C\\ is constructed to form
\\L-1\\ linearly independent comparisons among \\L\\ regimes, and the
chi-squared statistic is computed as \\(CZ)' (C \Sigma C')^{-1} (CZ)\\.
