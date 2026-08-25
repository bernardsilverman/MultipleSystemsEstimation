# Preliminary processing of a data matrix

Converts capture history data to the encoded representation used by the
model-fitting and extended-MLE routines.

## Usage

``` r
ingest_data(xdat)
```

## Arguments

- xdat:

  Capture history data in the standard package format.

## Value

A list with components:

- `nobs`:

  Counts indexed by encoded capture history.

- `nstar`:

  For each encoded history, the total count for that history and all its
  descendants.

- `nlists`:

  Number of capture lists.

- `listnames`:

  List names, constructed as A, B, and so on if necessary.

- `data`:

  The input data.

- `notestimable`:

  Logical vector identifying parameters that are strict descendants of
  parameters having zero sufficient statistic.

- `masterdesign`:

  The inclusion matrix constructed by
  [`make_master_design()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/make_master_design.md).

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection. *Statistics
and Computing*, **34**, 44.
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9)
.
