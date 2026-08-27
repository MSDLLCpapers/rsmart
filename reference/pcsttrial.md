# pcsttrial: Simulated Clinical Trial Data for Pain Coping Skills Training

A dataset containing simulated patient data to mimic the characteristics
of the Pain Coping Skills Training trial presented in Manschot, Laber,
and Davidian (2023).

## Usage

``` r
pcsttrial
```

## Format

A data frame with 284 rows and 14 variables:

- pctchange:

  Percent change in outcome from baseline at final assessment

- height:

  Patient height (cm)

- weight:

  Patient weight (kg)

- comorbidity:

  Comorbidity indicator (0 = no, 1 = yes)

- painmed:

  Pain medication use indicator (0 = no, 1 = yes)

- chemo:

  Chemotherapy indicator (0 = no, 1 = yes)

- pctchangek2:

  Percent change in outcome at stage 2

- adherence:

  Treatment adherence measure, from 0.5 to 1 in increments of 0.1

- a1:

  First-stage treatment assignment

- a2:

  Second-stage treatment assignment

- r2:

  Response indicator at stage 2 (0 = non-responder, 1 = responder)

- study_day_enroll:

  Study day of enrollment

- study_day_rerand:

  Study day of re-randomization

- study_day_outcome:

  Study day of outcome assessment

## Details

Note: this dataset is for illustrative purposes and does not use the
actual trial data.
