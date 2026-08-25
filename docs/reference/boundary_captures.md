# Find boundary terms of a hierarchical model

Finds encoded terms that are absent from the current hierarchical model
but whose immediate parents are all present. Adding any returned term
therefore preserves hierarchy.

## Usage

``` r
boundary_captures(kcap, nlists)
```

## Arguments

- kcap:

  Numeric vector containing the complete encoded parameter set of a
  hierarchical model.

- nlists:

  Total number of lists.

## Value

A numeric vector containing the admissible encoded boundary terms.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34**, 44,
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9)
.
