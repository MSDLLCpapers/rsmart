# Generate regime lists for a two-stage SMART where responders are not re-randomized

Constructs the treatment assignment matrices and consistency indicator
matrices for each embedded regime in a two-stage SMART. For stages where
response status determines treatment deterministically, the regime
matrix is updated to reflect the responder treatment.

## Usage

``` r
regime_list_no_trt_resp(
  emb_regimes,
  dat,
  resp_trt = list(r2 = list(0, 0, 0, 0))
)
```

## Arguments

- emb_regimes:

  A list of numeric vectors, each of length K, specifying the treatment
  codes (0 or 1) for non-responders at each stage. For example,
  `list(c(0,0), c(0,1), c(1,0), c(1,1))` encodes four regimes.

- dat:

  A data frame containing the SMART data. Expected to contain columns
  `a1, ..., aK` (treatment assignments) and optionally `r1, ..., rK`
  (response indicators). Typically this is the output from
  [`gen_no_trt_resp()`](https://msdllcpapers.github.io/rsmart/reference/gen_no_trt_resp.md),
  but any data frame with the required columns may be used.

- resp_trt:

  A named list where each element corresponds to a stage with
  deterministic responder treatment (e.g., `"r2"`). Each element is a
  list of length equal to the number of regimes, specifying the
  treatment assignment (0 or 1) for responders at that stage under each
  regime. Default is `list("r2" = list(0, 0, 0, 0))`.

## Value

A list of length `length(emb_regimes)`. Each element is a list with two
components:

- regime:

  An n x K matrix of treatment assignments each individual would receive
  if they followed that regime.

- regime_ind:

  An n x K indicator matrix (0 or 1) of whether each individual's
  observed treatment was consistent with the regime at each stage.
