# Simulate treatment assignments for a single stage

Generates treatment assignments for `n` individuals given a specified
number of possible treatments and a randomization probability function.
This is a utility function for constructing SMART simulations with
flexible randomization schemes. If an input data frame is provided, the
treatment assignments are appended as a new column `a<stage>`.

## Usage

``` r
sim_treatment(
  n = NULL,
  n_treatments,
  rand_prob_fn = function(n, n_treatments, prob) {
     sample(0:(n_treatments - 1), size
    = n, replace = TRUE, prob = prob)
 },
  dat = NULL,
  stage = 1,
  randomize_response = NULL,
  prob = NULL
)
```

## Arguments

- n:

  A positive integer. Number of individuals to assign treatments to. If
  `dat` is provided, `n` is inferred from `nrow(dat)` and this argument
  is ignored.

- n_treatments:

  A positive integer. Number of possible treatments. Treatments are
  coded as integers `0, 1, ..., n_treatments - 1`.

- rand_prob_fn:

  A function that takes `n` (the number of individuals), `n_treatments`
  (the number of treatment options), and `prob` (a numeric vector of
  probabilities) as arguments, and returns an integer vector of length
  `n` with treatment assignments coded as `0, 1, ..., n_treatments - 1`.
  Default is Bernoulli randomization using the probabilities specified
  by `prob`.

- dat:

  An optional data frame. If provided, the function appends the
  treatment assignments as a new column named `a<stage>` and returns the
  modified data frame. The data frame must not already contain a column
  with that name. Required when `stage > 1`.

- stage:

  A positive integer indicating the stage number. The treatment column
  will be named `paste0("a", stage)`. Default is 1. When `stage > 1`,
  `dat` must be provided and must contain columns `a1, ..., a<stage-1>`.
  Randomization is performed independently within each unique
  combination of prior treatments.

- randomize_response:

  One of `NULL` (default), `"Y"`, or `"N"`. If `"Y"`, the response
  indicator `r<stage>` is included in the grouping variable along with
  prior treatments `a1, ..., a<stage-1>` when randomizing. This allows
  different randomization within responder and non-responder subgroups.
  If `"N"`, responders (`r<stage> == 1`) are assigned treatment `0`
  deterministically, and only non-responders are randomized within prior
  treatment groups. Both `"Y"` and `"N"` require `stage > 1` and `dat`
  to contain a column named `r<stage>`. If `NULL`, response status is
  not used.

- prob:

  A numeric vector or a named list specifying randomization
  probabilities. If a numeric vector, it is used as the probability
  weights for all groups (must have length equal to `n_treatments`). If
  a named list, the names should correspond to the group levels (formed
  by [`interaction()`](https://rdrr.io/r/base/interaction.html) of prior
  treatment columns and optionally the response column), and each
  element should be a numeric vector of length `n_treatments` giving the
  group-specific probabilities. Default is `NULL`, which uses equal
  probability across treatments.

## Value

If `dat` is `NULL` (default), an integer vector of length `n` with
treatment assignments coded as `0, 1, ..., n_treatments - 1`. If `dat`
is provided, the input data frame with a new integer column `a<stage>`
appended.

## Examples

``` r
# Default: Bernoulli randomization with 2 treatments (prob 0.5)
set.seed(1)
a <- sim_treatment(n = 100, n_treatments = 2)
table(a)
#> a
#>  0  1 
#> 48 52 

# Unequal randomization: 70% to treatment 0, 30% to treatment 1
set.seed(1)
a <- sim_treatment(n = 100, n_treatments = 2, prob = c(0.7, 0.3))
table(a)
#> a
#>  0  1 
#> 68 32 

# Three treatments with equal probability
set.seed(1)
a <- sim_treatment(n = 300, n_treatments = 3)
table(a)
#> a
#>   0   1   2 
#>  90 101 109 

# With an input data frame
set.seed(1)
df <- data.frame(x1 = rnorm(100), x2 = rbinom(100, 1, 0.5))
df <- sim_treatment(n_treatments = 2, dat = df)
head(df)
#>           x1 x2 a1
#> 1 -0.6264538  0  0
#> 2  0.1836433  0  1
#> 3 -0.8356286  1  1
#> 4  1.5952808  0  1
#> 5  0.3295078  0  1
#> 6 -0.8204684  1  0

# With an input data frame at stage 2
set.seed(1)
df <- data.frame(x1 = rnorm(100), a1 = rbinom(100, 1, 0.5))
df <- sim_treatment(n_treatments = 2, dat = df, stage = 2)
head(df)
#>           x1 a1 a2
#> 1 -0.6264538  0  0
#> 2  0.1836433  0  1
#> 3 -0.8356286  1  0
#> 4  1.5952808  0  1
#> 5  0.3295078  0  1
#> 6 -0.8204684  1  1

if (FALSE) { # \dontrun{
# Stage 2 with different randomization probabilities depending on a1:
# equal (50/50) if a1 == 0, unequal (70/30) if a1 == 1
# Using a named list for prob, where names match group levels
set.seed(1)
n <- 200
df <- data.frame(a1 = rbinom(n, 1, 0.5))
df <- sim_treatment(n_treatments = 2, dat = df, stage = 2,
                    prob = list("0" = c(0.5, 0.5), "1" = c(0.7, 0.3)))
table(df$a1, df$a2)
} # }
```
