# Jackknife abundance and Jackknife bic

This routine takes the output from `subsetmat` or from `assemble_bic`
and returns the jackknife abundance matrix and jackknife BIC matrix.

## Usage

``` r
jackknifecal(z, checkexist = TRUE)
```

## Arguments

- z:

  Results from `assemble_bic` or `subsetmat`.

- checkexist:

  If `checkexist=TRUE`, check for existence in cases where the jackknife
  introduces an additional zero, else it does not check for existence.
  Note that in the current version it is assume that models for which
  the fit doesn't exist for the original data have already been
  excluded.

## Value

A list with the following components

- jackabund:

  Jackknife abundance matrix

- jackbic:

  Jackknife BIC matrix

- countsobserved:

  Capture counts in the same order as the columns of `jackabund` and
  `jackbic`

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
[\doi{10.1007/s11222-023-10346-9}](NA).
