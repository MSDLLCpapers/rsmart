# Fit an outcome regression (Q-function) model for a single stage

Fits a single Q-function model at one stage of the SMART. Predictions
are obtained both for the unmodified data and for data where the
treatment is set to the recommended regime.

## Usage

``` r
qstep(qmodel, data, response, newdata, regime, txName)
```

## Arguments

- qmodel:

  A `modelObj` object specifying the outcome regression model.

- data:

  A data frame of individuals used to fit the model (typically those who
  have completed the stage).

- response:

  A numeric vector or single-column data frame of responses
  (pseudo-outcomes from later stages or observed outcomes).

- newdata:

  A data frame of individuals for whom predictions are desired
  (typically those who have reached the stage).

- regime:

  A vector of recommended treatment assignments under the regime being
  evaluated.

- txName:

  A character string specifying the name of the treatment column to
  modify (e.g., `"a1"`, `"a2"`).

## Value

A list with the following components:

- hats_mod:

  A numeric vector of predicted values under the regime-consistent
  treatment.

- hats_unmod:

  A numeric vector of predicted values under the actual (unmodified)
  treatment.

- qfit:

  The fitted `modelObj` object.
