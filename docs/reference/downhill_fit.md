# Downhill search among hierarchical models

Starting from the main-effects model, repeatedly moves to the
neighbouring hierarchical model with the smallest BIC until no
improvement is available or the iteration limit is reached.

## Usage

``` r
downhill_fit(
  counts,
  desmat,
  maxorder = dim(desmat)[2] - 1,
  checkid = TRUE,
  niter = 20,
  verbose = FALSE
)
```

## Arguments

- counts:

  Numeric vector of observed counts.

- desmat:

  Binary matrix defining the corresponding capture histories.

- maxorder:

  Maximum interaction order considered.

- checkid:

  If `TRUE`, check parameter identifiability and existence of the
  extended MLE before fitting each model.

- niter:

  Maximum number of downhill iterations.

- verbose:

  If `FALSE`, return only the selected population estimate. If `TRUE`,
  return detailed search results.

## Value

If `verbose = FALSE`, the population estimate from the selected model,
or `NA` if the main-effects model has no valid fit. If `verbose = TRUE`,
a list with components:

- `optimum_hierarchy`:

  Selected hierarchical model.

- `minimum_value`:

  Named vector containing its BIC and population estimate.

- `hierarchies_considered`:

  Character vector of models examined.

- `function_values`:

  Matrix containing the BIC and population estimate for each model
  examined.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34**, 44.
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9)
.
