# Jackknife abundance and BIC values

Constructs the delete-one jackknife fits needed for BCa acceleration.
Each distinct positive-count capture history is reduced by one
individual in turn, and every retained hierarchical model is fitted to
the resulting data.

## Usage

``` r
jackknifecal(z, checkexist = TRUE)
```

## Arguments

- z:

  A result from
  [`assemble_bic()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/assemble_bic.md)
  or
  [`subsetmat()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/subsetmat.md).

- checkexist:

  If `TRUE`, check identifiability and existence of the extended MLE
  when a deletion creates an additional zero count. Models that fail on
  the original data are assumed to have been removed already.

## Value

The input list `z`, with the following components added or replaced:

- `jackabund`:

  A matrix of jackknife population estimates, with models in rows and
  capture histories in columns.

- `jackbic`:

  A corresponding matrix of BIC values.

- `countsobserved`:

  Capture counts in the same order as the columns of `jackabund` and
  `jackbic`.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection. *Statistics
and Computing*, **34**, 44.
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9).
