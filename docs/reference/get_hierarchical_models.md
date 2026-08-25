# Extract hierarchical models from a catalogue

Selects all hierarchical models for a specified number of lists and
maximum interaction order from a catalogue of model strings.

## Usage

``` r
get_hierarchical_models(nlists, maxorder = nlists - 1, modelvec = hiermodels)
```

## Arguments

- nlists:

  Number of lists.

- maxorder:

  Maximum interaction order. The default is `nlists - 1`.

- modelvec:

  A character vector containing the catalogue of hierarchical-model
  strings from which models are selected. The default is the precomputed
  package catalogue `hiermodels`. An alternative catalogue using the
  same hierarchy-string notation may be supplied. The default catalogue
  contains all hierarchical models for two to five lists, and all
  six-list hierarchical models with interaction order at most 2.

## Value

A character vector containing the models satisfying the specified
number-of-lists and maximum-interaction-order restrictions.

## References

Silverman, B. W., Chan, L. and Vincent, K. (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection. *Statistics
and Computing*, **34**, 44.
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9)
.

## Examples

``` r
# Three lists, all interaction orders
get_hierarchical_models(nlists = 3, maxorder = Inf)
#> [1] "[12,13,23]" "[12,13]"    "[12,23]"    "[13,23]"    "[12,3]"    
#> [6] "[13,2]"     "[23,1]"     "[1,2,3]"   

# Three lists with the default maximum interaction order of 2
get_hierarchical_models(nlists = 3)
#> [1] "[12,13,23]" "[12,13]"    "[12,23]"    "[13,23]"    "[12,3]"    
#> [6] "[13,2]"     "[23,1]"     "[1,2,3]"   

# Number of five-list models with maximum interaction order 3
length(get_hierarchical_models(nlists = 5, maxorder = 3))
#> [1] 6212
```
