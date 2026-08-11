# Find the vector of captures corresponding to a given hierarchical model

Given a hierarchical model, find the vector of all the corresponding
encoded captures

## Usage

``` r
convert_from_hierarchy(modelstr, findancestors = TRUE)
```

## Arguments

- modelstr:

  A given hierarchical model

- findancestors:

  If TRUE then find all the captures. If FALSE then just return the
  encoded defining histories of the hierarchy

## Value

The encoded capture histories that corresponds to the row number of the
capture history data set

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
[\doi{10.1007/s11222-023-10346-9}](NA).

## Examples

``` r
modelstr = "[12,23]"
convert_from_hierarchy(modelstr)
#> [1] 1 2 3 4 5 7
modelstr = "[12,3]"
convert_from_hierarchy(modelstr, findancestors=FALSE)
#> [1] 4 5
```
