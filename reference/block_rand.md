# Create a blocked randomization function

Returns a randomization function compatible with
[`sim_treatment`](https://msdllcpapers.github.io/rsmart/reference/sim_treatment.md)
that implements permuted block randomization. Within each block, every
treatment appears exactly `block_rep` times, and the order within blocks
is randomly permuted. This ensures approximate balance across treatments
at any point during enrollment.

## Usage

``` r
block_rand(block_rep = 2)
```

## Arguments

- block_rep:

  A positive integer. Number of times each treatment appears in every
  block. The block size is `block_rep * n_treatments`. Default is 2.

## Value

A function with signature `function(n, n_treatments, prob)` suitable for
use as `rand_prob_fn` in
[`sim_treatment`](https://msdllcpapers.github.io/rsmart/reference/sim_treatment.md).
The `prob` argument is accepted but ignored, since blocked randomization
enforces equal allocation within each block.

## Examples

``` r
# Use blocked randomization with sim_treatment
set.seed(1)
a <- sim_treatment(n = 100, n_treatments = 2,
                   rand_prob_fn = block_rand(block_rep = 2))
table(a)
#> a
#>  0  1 
#> 50 50 

# Three treatments with block size 6 (block_rep = 2)
set.seed(1)
a <- sim_treatment(n = 300, n_treatments = 3,
                   rand_prob_fn = block_rand(block_rep = 2))
table(a)
#> a
#>   0   1   2 
#> 100 100 100 

# Stage 2 blocked randomization within prior treatment groups
set.seed(1)
df <- data.frame(a1 = c(rep(0, 50), rep(1, 50)))
df <- sim_treatment(n_treatments = 2, dat = df, stage = 2,
                    rand_prob_fn = block_rand(block_rep = 3))
table(df$a1, df$a2)
#>    
#>      0  1
#>   0 24 26
#>   1 26 24
```
