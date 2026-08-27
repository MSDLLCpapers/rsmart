# pkg_architecture

``` r

library(rsmart)
```

The goal of this vignette is to help users understand the overall
structure of the **rsmart** package. While the ‘Reference’ page on the
pkgdown website is helpful for examining functions individually, we want
to provide a high-level view of the package’s conceptual structure.

The **rsmart** package has **3 user-facing functions**, along with a
hierarchy of internal helpers. These functions can be organized into
five conceptual groups, each of which we will unpack here in this
vignette.

### 1. User-Facing Entry Points

- **[`iaipwe()`](https://msdllcpapers.github.io/rsmart/reference/iaipwe.md)**
  — The main estimation workhorse. Orchestrates the entire pipeline:
  computes kappa, nu, propensity scores, regime values, and then the
  sandwich variance (or bootstrap variance). This is what users call
  directly.
- **[`gen_no_trt_resp()`](https://msdllcpapers.github.io/rsmart/reference/gen_no_trt_resp.md)**
  — Data generation for simulation studies.
- **[`regime_list_no_trt_resp()`](https://msdllcpapers.github.io/rsmart/reference/regime_list_no_trt_resp.md)**
  — Builds the regime/indicator matrices that
  [`iaipwe()`](https://msdllcpapers.github.io/rsmart/reference/iaipwe.md)
  needs.

### 2. Nuisance Parameter Estimation (called by `iaipwe`)

These are called *first* inside
[`iaipwe()`](https://msdllcpapers.github.io/rsmart/reference/iaipwe.md)
to estimate the building blocks:

- **[`get_kappa()`](https://msdllcpapers.github.io/rsmart/reference/get_kappa.md)**
  — Computes how far each individual progressed (stage reached).
- **[`get_nu()`](https://msdllcpapers.github.io/rsmart/reference/get_nu.md)**
  — Estimates stage-arrival probabilities \\\nu_k\\.
- **[`pi_fits()`](https://msdllcpapers.github.io/rsmart/reference/pi_fits.md)**
  — Fits propensity models at all stages (calls
  **[`pstep()`](https://msdllcpapers.github.io/rsmart/reference/pstep.md)**
  per stage).

### 3. Value Estimation (called by `iaipwe`)

- **[`estimate_values()`](https://msdllcpapers.github.io/rsmart/reference/estimate_values.md)**
  — Loops over regimes and for each one:
  - **[`get_q_fits()`](https://msdllcpapers.github.io/rsmart/reference/get_q_fits.md)**
    — Fits Q-functions backwards through stages (calls
    **[`qstep()`](https://msdllcpapers.github.io/rsmart/reference/qstep.md)**
    per stage).
  - **[`value_terms()`](https://msdllcpapers.github.io/rsmart/reference/value_terms.md)**
    — Computes the \\2K+1\\ coarsening-level value terms (augmentation +
    IPW).

### 4. Sandwich Variance — \\B_n\\ (empirical variance of estimating equations)

**[`get_bn()`](https://msdllcpapers.github.io/rsmart/reference/get_bn.md)**
assembles \\\Psi_i \Psi_i^T / n\\ by collecting individual-level
estimating equation contributions:

- **`ee_psi_pi()`** — Estimating equation contributions for \\\pi\\
  parameters.
- **`ee_psi_nu()`** — Estimating equation contributions for \\\nu\\
  parameters.
- **`ee_psi_beta()`** — Estimating equation contributions for \\\beta\\
  (Q-function) parameters.
- **`ee_psi_v()`** — Estimating equation contributions for the value
  parameters \\V\\.

### 5. Sandwich Variance — \\A_n\\ (derivative of estimating equations)

**[`get_an()`](https://msdllcpapers.github.io/rsmart/reference/get_an.md)**
assembles \\-\partial\Psi/\partial\theta\\ by collecting derivatives:

- **`ee_dpsi_pi()`** — Derivative block for \\\pi\\.
- **`ee_dpsi_nu()`** — Derivative block for \\\nu\\.
- **`ee_dpsi_beta()`** — Derivative block for \\\beta\\.
- **`ee_dpsiv()`** (AIPW) or **`ee_dpsiv_ipw()`** (IPW) — Derivative
  rows for \\V\\, which internally call:
  - **`ee_dpsiv_dpi()`** / **`ee_dpsiv_dpi_ipw()`**
  - **`ee_dpsiv_dnu()`** / **`ee_dpsiv_dnu_ipw()`**
  - **`ee_dpsiv_dbeta()`** (AIPW only)
  - **`ee_dpsiv_dv()`**

### 6. Trial Design Utilities (standalone)

- **[`get_bounds()`](https://msdllcpapers.github.io/rsmart/reference/get_bounds.md)**
  — Group sequential stopping boundaries (calls
  **[`get_first_bound()`](https://msdllcpapers.github.io/rsmart/reference/get_first_bound.md)** +
  **[`get_next_bound()`](https://msdllcpapers.github.io/rsmart/reference/get_next_bound.md)**).
- **[`get_sample_size()`](https://msdllcpapers.github.io/rsmart/reference/get_sample_size.md)**
  — Sample size determination.
- **[`get_q_coefs()`](https://msdllcpapers.github.io/rsmart/reference/get_q_coefs.md)**
  — Coefficient extraction utility.
