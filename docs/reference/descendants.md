# Find descendants of an encoded capture history

Finds all encoded histories that contain the supplied history.

## Usage

``` r
descendants(k, nlists, omitk = FALSE)
```

## Arguments

- k:

  An encoded capture history.

- nlists:

  Total number of lists.

- omitk:

  If `TRUE`, omit `k` itself from the result.

## Value

A sorted numeric vector containing the encoded descendants.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34**, 44,
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9).
