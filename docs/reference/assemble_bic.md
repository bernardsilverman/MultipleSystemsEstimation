# Enumerate and rank hierarchical models by BIC

Fits the hierarchical models permitted by the specified number of lists
and maximum interaction order, and orders the valid fits by increasing
BIC.

## Usage

``` r
assemble_bic(
  xdata,
  maxorder = dim(xdata)[2] - 2,
  checkexist = TRUE,
  removeFRfail = TRUE,
  ...
)
```

## Arguments

- xdata:

  Capture history data in the standard package format.

- maxorder:

  Maximum interaction order to be included.

- checkexist:

  If `TRUE`, check parameter identifiability and existence of the
  extended MLE for each model.

- removeFRfail:

  If `TRUE`, remove models without a valid fit from the results.

- ...:

  Additional arguments passed to
  [`fit_hier_model()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/fit_hier_model.md).

## Value

A list with components:

- `res`:

  A matrix containing the population estimate, BIC and maximum
  interaction order for each retained model, ordered by increasing BIC.
  Model strings are used as row names.

- `xdata`:

  The original capture history data.

- `maxorder`:

  The largest interaction order among the retained models, or 0 if no
  model is retained.

- `best_neginfpars`:

  Encoded effects estimated at minus infinity in the model with the
  smallest BIC, or an empty vector if there is no valid model or no such
  effect.
