# Package index

## Primary Analysis

Main entry point for IAIPWE estimation

- [`iaipwe()`](https://msdllcpapers.github.io/rsmart/reference/iaipwe.md)
  : IAIPWE for K-stage SMARTs with up to 2 treatment options at each
  stage
- [`estimate_values()`](https://msdllcpapers.github.io/rsmart/reference/estimate_values.md)
  : Estimate values for all treatment regimes

## Trial Design

Functions for SMART design, sample size, and stopping boundaries

- [`smart_design()`](https://msdllcpapers.github.io/rsmart/reference/smart_design.md)
  : Compute stopping boundaries and sample size for a group sequential
  SMART
- [`get_sample_size()`](https://msdllcpapers.github.io/rsmart/reference/get_sample_size.md)
  : Determine sample size for a group sequential SMART design
- [`get_sample_size_chi()`](https://msdllcpapers.github.io/rsmart/reference/get_sample_size_chi.md)
  : Determine sample size for a chi-squared global test in a group
  sequential SMART design
- [`get_bounds()`](https://msdllcpapers.github.io/rsmart/reference/get_bounds.md)
  : Compute all stopping boundaries for a two-analysis group sequential
  design
- [`get_bounds_chi()`](https://msdllcpapers.github.io/rsmart/reference/get_bounds_chi.md)
  : Compute the first chi-squared stopping boundary for a group
  sequential design
- [`get_first_bound()`](https://msdllcpapers.github.io/rsmart/reference/get_first_bound.md)
  : Compute the first stopping boundary for a group sequential design
- [`get_next_bound()`](https://msdllcpapers.github.io/rsmart/reference/get_next_bound.md)
  : Compute subsequent stopping boundaries for a group sequential design

## Data Generation and Simulation

Helpers for generating treatment regimes and simulating trials

- [`gen_no_trt_resp()`](https://msdllcpapers.github.io/rsmart/reference/gen_no_trt_resp.md)
  : Generate sample data for a two-stage SMART where responders are not
  re-randomized
- [`regime_list_no_trt_resp()`](https://msdllcpapers.github.io/rsmart/reference/regime_list_no_trt_resp.md)
  : Generate regime lists for a two-stage SMART where responders are not
  re-randomized
- [`sim_treatment()`](https://msdllcpapers.github.io/rsmart/reference/sim_treatment.md)
  : Simulate treatment assignments for a single stage
- [`block_rand()`](https://msdllcpapers.github.io/rsmart/reference/block_rand.md)
  : Create a blocked randomization function

## Model Fitting

Q-learning and propensity score model fitting

- [`qstep()`](https://msdllcpapers.github.io/rsmart/reference/qstep.md)
  : Fit an outcome regression (Q-function) model for a single stage
- [`pstep()`](https://msdllcpapers.github.io/rsmart/reference/pstep.md)
  : Fit a propensity score model for a single stage
- [`get_q_fits()`](https://msdllcpapers.github.io/rsmart/reference/get_q_fits.md)
  : Fit outcome regression (Q-function) models across all stages for a
  single regime
- [`get_q_coefs()`](https://msdllcpapers.github.io/rsmart/reference/get_q_coefs.md)
  : Extract coefficients from all fitted Q-function models
- [`pi_fits()`](https://msdllcpapers.github.io/rsmart/reference/pi_fits.md)
  : Fit propensity score models for all stages

## Variance Estimation

Sandwich variance components

- [`get_an()`](https://msdllcpapers.github.io/rsmart/reference/get_an.md)
  : Compute the An matrix for the sandwich variance estimator
- [`get_bn()`](https://msdllcpapers.github.io/rsmart/reference/get_bn.md)
  : Compute the Bn matrix for the sandwich variance estimator
- [`value_terms()`](https://msdllcpapers.github.io/rsmart/reference/value_terms.md)
  : Compute the individual-level value terms for a single regime
- [`get_nu()`](https://msdllcpapers.github.io/rsmart/reference/get_nu.md)
  : Estimate stage arrival probabilities (nu)
- [`get_kappa()`](https://msdllcpapers.github.io/rsmart/reference/get_kappa.md)
  : Compute the kappa (stage reached) for each individual

## Data

Example datasets

- [`pcsttrial`](https://msdllcpapers.github.io/rsmart/reference/pcsttrial.md)
  : pcsttrial: Simulated Clinical Trial Data for Pain Coping Skills
  Training
