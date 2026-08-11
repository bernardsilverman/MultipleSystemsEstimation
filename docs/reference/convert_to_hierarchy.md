# Find hierarchical representation of a vector of captures

Given a vector of encoded captures defining a hierarchical model,
re-express it in hierarchical model form

## Usage

``` r
convert_to_hierarchy(kcap, nlists)
```

## Arguments

- kcap:

  A vector of captures

- nlists:

  The number of lists

## Value

A hierarchical representation of the vector of encoded captures.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
<https://doi.org/10.1007/s11222-023-10346-9>.

## Examples

``` r
kcap=c(1,2,3,5,4)
nlists=3
convert_to_hierarchy(kcap, nlists)
#> [1] "[12,3]"
```
