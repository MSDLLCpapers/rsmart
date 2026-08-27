# Fit a propensity score model for a single stage

Fits a propensity score model at a single stage of the SMART and returns
the fitted probabilities. The probability that each individual received
the treatment they actually received is computed.

## Usage

``` r
pstep(pmodel, data, response, k)
```

## Arguments

- pmodel:

  A `modelObj` object specifying the propensity score model (typically a
  binomial GLM).

- data:

  A data frame of individuals to fit the model on (typically those who
  have reached stage `k`).

- response:

  A vector or data frame column of treatment assignments at stage `k`.

- k:

  An integer indicating the stage number.

## Value

A list with the following components:

- pk:

  A numeric vector of predicted probabilities that treatment is 1 at
  stage `k`.

- ps:

  A numeric vector of estimated propensity scores (probability that the
  individual received the treatment they actually received).

- pfit:

  The fitted `modelObj` object.
