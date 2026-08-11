# Find the "ancestors" of a given capture history

Given any encoded capture history and the number of lists, find all the
encoded capture histories that are included in the original capture
history

## Usage

``` r
ancestors(k, nlists = 10)
```

## Arguments

- k:

  An encoded capture history

- nlists:

  The total number of lists

## Value

a vector giving the encoded versions of the ancestors

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
<https://doi.org/10.1007/s11222-023-10346-9>.

## Examples

``` r
ancestors(2,10)
#> [1] 1 2
ancestors(1,5)
#> [1] 1
```
