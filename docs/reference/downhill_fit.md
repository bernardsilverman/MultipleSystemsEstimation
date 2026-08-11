# Conduct downhill search among hierarchical models starting from the main effects only.

Find a local optimum by downhill search among hierarchical models

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

  Observed counts for the capture histories defined by desmat

- desmat:

  Incidence matrix defining the capture histories observed with counts
  given by counts

- maxorder:

  Maximum order of models to be included

- checkid:

  If it is TRUE, then `checkident.1` is called and it performs the
  Fienberg-Rinaldo linear program check for the existence of the
  estimates

- niter:

  Number of iterations

- verbose:

  Specifies the output, if FALSE then only returns the best value, if
  TRUE, returns a more detailed list of objects

## Value

A list with the following components

- optimum_hierarchy:

  Optimal hierarchical model

- minimum_value:

  hierarchical model with the minimum value

- hierarchies_considered:

  hierarhical models considered

- function_values:

  Values of function

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
<https://doi.org/10.1007/s11222-023-10346-9>.
