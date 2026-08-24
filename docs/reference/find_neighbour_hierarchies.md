# Find the neighbours of a hierarchical model

Given a hierarchical model, finds its outer neighbours, its inner
neighbours, or both.

## Usage

``` r
find_neighbour_hierarchies(
  modelstr,
  nlists = NA,
  keepmaineffects = TRUE,
  maxorder = nlists - 1,
  type = c("all", "outer", "inner")
)
```

## Arguments

- modelstr:

  A model string written in hierarchical form.

- nlists:

  The total number of lists. If `NA`, it is inferred from the largest
  list number appearing in `modelstr`.

- keepmaineffects:

  If `TRUE`, main effects are not removed when constructing inner
  neighbours.

- maxorder:

  The maximum interaction order permitted in outer neighbours.

- type:

  Which neighbours to return. `"all"` returns both inner and outer
  neighbours, `"outer"` returns only models obtained by adding a term,
  and `"inner"` returns only models obtained by removing a generator.

## Value

A character vector containing the requested neighbouring hierarchical
models.

## Details

An outer neighbour is obtained by adding an interaction term all of
whose subsets are already in the model.

An inner neighbour is obtained by removing one of the generators
defining the hierarchical model. Removing a generator removes only the
defining interaction term itself, not the lower-order terms that it
implies. For example, removing the generator `123` from `[123,34]`
yields `[12,13,23,34]`, not `[34]`. Main effects are not removed when
`keepmaineffects = TRUE`.

## References

Silverman, B. W., Chan, L. and Vincent, K. (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection. *Statistics
and Computing*, **34**, 44.
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9).

## Examples

``` r
modelstr <- "[12,23]"

find_neighbour_hierarchies(modelstr)
#> [1] "[23,1]"     "[12,3]"     "[12,13,23]"
find_neighbour_hierarchies(modelstr, type = "outer")
#> [1] "[12,13,23]"
find_neighbour_hierarchies(modelstr, type = "inner")
#> [1] "[23,1]" "[12,3]"
```
