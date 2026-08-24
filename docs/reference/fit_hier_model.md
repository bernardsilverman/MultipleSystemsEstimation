# Fit a hierarchical model taking account of possible sparsity

Fits a Poisson log-linear model on the appropriate face of the parameter
space. Parameters whose sufficient statistics are zero are fixed at
minus infinity and their descendant cells are removed before fitting.

## Usage

``` r
fit_hier_model(xdatin, hiermod, bicRcap = TRUE, checkid = FALSE)
```

## Arguments

- xdatin:

  Data prepared by
  [`ingest_data()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/ingest_data.md).

- hiermod:

  Character string specifying the hierarchical model.

- bicRcap:

  If `TRUE`, use the number of observed cases as the BIC sample size.
  Otherwise use the number of fitted cells in the Poisson log-linear
  model.

- checkid:

  If `TRUE`, check parameter identifiability and existence of the
  extended MLE before fitting.

## Value

An object returned by
[`stats::glm.fit()`](https://rdrr.io/r/stats/glm.html), augmented by:

- `abundance`:

  Estimated total population size, or `NA` if the model has no valid
  fit.

- `bic`, `aic`:

  BIC and AIC values, where available.

- `neginfpars`:

  Encoded parameters whose extended-MLE values are minus infinity.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection. *Statistics
and Computing*, **34**, 44.
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9).
