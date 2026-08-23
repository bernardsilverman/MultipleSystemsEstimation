# Find the "descendants" of a given capture history

Given any encoded capture history, find all the encoded capture
histories that include the original capture history and any other lists

## Usage

``` r
descendants(k, nlists, omitk = FALSE)
```

## Arguments

- k:

  An encoded capture history

- nlists:

  The total number of lists

- omitk:

  Determine whether the original capture history is included as a
  descendant of itself. If `omitk=TRUE` it is not.

## Value

a vector giving the encoded versions of the descendants

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
[\doi{10.1007/s11222-023-10346-9}](NA).

## Examples

``` r
descendants(2,5)
#> Error in descendants(2, 5): could not find function "descendants"
descendants(5,10)
#> Error in descendants(5, 10): could not find function "descendants"
```
