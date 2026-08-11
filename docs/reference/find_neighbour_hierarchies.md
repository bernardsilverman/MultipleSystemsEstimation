# Find all neighbouring hierarchical model to a given one

Find all neighbouring hierarchical model to a given one

## Usage

``` r
find_neighbour_hierarchies(
  modelstr,
  nlists = NA,
  keepmaineffects = TRUE,
  maxorder = nlists - 1
)
```

## Arguments

- modelstr:

  Model string written in hierarchical form

- nlists:

  Number of lists.

- keepmaineffects:

  If TRUE, keep the main effects. If FALSE remove.

- maxorder:

  Maximum order of models to be included

## Value

neighbour hierarchical models

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
[\doi{10.1007/s11222-023-10346-9}](NA).

## Examples

``` r
modelstr = "[12,23]"
find_neighbour_hierarchies(modelstr)
#> [1] "[23,1]"     "[12,3]"     "[12,13,23]"
```
