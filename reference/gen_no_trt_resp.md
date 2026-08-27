# Generate sample data for a two-stage SMART where responders are not re-randomized

This function generates data for a two-stage Sequential Multiple
Assignment Randomized Trial (SMART) where responders are not
re-randomized in the second stage. The function allows for different
value patterns and treatment assignments.

## Usage

``` r
gen_no_trt_resp(n, s2 = 100, block_rep = 2, r2p = 0.5)
```

## Arguments

- n:

  A positive integer. Number of individuals to generate data for.

- s2:

  A positive numeric value. Variance of the error term. Default is 100.

- block_rep:

  A positive integer. Number of times to duplicate treatments per block
  in the permuted block randomization. Default is 2.

- r2p:

  A numeric value between 0 and 1. Probability of being a responder in
  the second stage. Default is 0.5.

## Value

A data frame with `n` rows and the following columns:

- t1:

  Study day of enrollment.

- t2:

  Study day of second-stage randomization.

- t3:

  Study day of outcome assessment.

- a1:

  First-stage treatment assignment (0 or 1).

- r2:

  Response indicator at stage 2 (0 = non-responder, 1 = responder).

- a2:

  Second-stage treatment assignment (0 or 1 for non-responders, 0 for
  responders).

- x11:

  Baseline covariate (continuous, range 25–75).

- x12:

  Baseline covariate (binary, 0 or 1).

- x21:

  Stage 2 covariate (continuous, range 0–1).

- y:

  Continuous outcome (higher is better).

- id:

  Individual identifier (1 to `n`).

## Examples

``` r
set.seed(1)
dat <- gen_no_trt_resp(n=400, s2=100, block_rep=2, r2p = 0.5)
head(dat)
#>     t1  t2  t3 a1 r2 a2 x11 x12  x21        y id
#> 116 14 108 207  0  0  1  40   0 0.97 24.50311  1
#> 27  15 107 201  0  1  0  36   0 0.00 41.99025  2
#> 47  26 123 214  1  0  0  70   0 0.11 28.47873  3
#> 281 30 135 232  1  0  0  72   0 0.55 39.39817  4
#> 133 39 136 231  1  1  0  50   0 0.87 49.68450  5
#> 228 51 154 249  0  0  1  33   0 0.86 34.43490  6
```
