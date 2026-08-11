# Find the "parents" of a given capture history

Given any encoded capture history and the number of lists, find the
encoded capture histories which are obtained by leaving out just one
list in turn

## Usage

``` r
parent_captures(k, nlists = 10)
```

## Arguments

- k:

  An encoded capture history that corresponds to the row number of the
  capture history data set

- nlists:

  The total number of lists

## Value

a vector giving the encoded versions of the parents

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
<https://doi.org/10.1007/s11222-023-10346-9>.
