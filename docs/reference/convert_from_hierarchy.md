# Convert a hierarchy string to encoded parameters

Converts a hierarchical-model specification to its encoded generators or
to the complete hierarchical closure.

## Usage

``` r
convert_from_hierarchy(modelstr, findancestors = TRUE)
```

## Arguments

- modelstr:

  Character string specifying a hierarchical model.

- findancestors:

  If `TRUE`, return the complete hierarchical closure. If `FALSE`,
  return only the encoded generators.

## Value

A numeric vector containing the requested encoded parameters.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34**, 44,
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9).
