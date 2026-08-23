# Find the "ancestors" of a given capture history

Given any encoded capture history, find all the encoded capture
histories that are included in the original capture history

## Usage

``` r
ancestors(k)
```

## Arguments

- k:

  An encoded capture history

## Value

a vector giving the encoded versions of the ancestors

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
[\doi{10.1007/s11222-023-10346-9}](NA).

## Examples

``` r
ancestors(2)
#> Error in ancestors(2): could not find function "ancestors"
ancestors(1)
#> Error in ancestors(1): could not find function "ancestors"
```
