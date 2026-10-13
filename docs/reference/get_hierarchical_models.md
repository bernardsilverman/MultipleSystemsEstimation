# Generate hierarchical models

Generates all hierarchical models satisfying specified criteria. The
user can specify the number of lists, the maximum order of models to be
considered, and (optionally) a data set for which any masked sparse
interactions will not be included in the models. Intercept and main
effects are always included.

## Usage

``` r
get_hierarchical_models(
  nlists,
  maxorder = nlists - 1L,
  xdata = NULL,
  maxmodel = 100000L
)
```

## Arguments

- nlists:

  Number of lists (an integer from 2 to 9).

- maxorder:

  Maximum interaction order; Inf permits all orders.

- xdata:

  Optional capture history data. If supplied, models including masked
  sparse interactions for this data set will not be included. Minimal
  sparse interactions remain eligible.

- maxmodel:

  Maximum number of models to be generated. If exceeded, the maximum
  order will be reduced, if necessary as far as 1, to bring the number
  of models below the specified limit.

## Value

A named list of integer vectors specifying the hierarchical models. The
name of each vector in the list is the hierarchy string representation
of the model, and the vector itself is the vector of encoded capture
histories included in the model. A capture history is encoded as the
integer \\1 + \sum\_{i \in S} 2^{i-1}\\, where \\S\\ is the set of
capture lists in the history. The returned list has two attributes:

- `requested_maxorder`:

  The value of `maxorder` supplied by the user, or its default value if
  omitted.

- `maxorder`:

  The actual maximum interaction order of the generated models. This is
  capped at `nlists` and may be reduced further to satisfy `maxmodel`.

## Examples

``` r
models <- get_hierarchical_models(3)
names(models)
#> [1] "[1,2,3]"    "[23,1]"     "[13,2]"     "[13,23]"    "[12,3]"    
#> [6] "[12,23]"    "[12,13]"    "[12,13,23]"
models[[1]]
#> [1] 1 2 3 5
```
