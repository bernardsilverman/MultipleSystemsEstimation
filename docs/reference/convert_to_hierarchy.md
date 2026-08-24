# Find hierarchical representation of a vector of captures

Given a vector of encoded captures defining a hierarchical model,
re-express it in hierarchical model form.

## Usage

``` r
convert_to_hierarchy(kcap)
```

## Arguments

- kcap:

  Numeric vector of encoded parameters.

## Value

A hierarchical representation of the vector of encoded captures.

## Details

The supplied parameter vector must define a hierarchical model: whenever
an interaction is present, all its lower-order terms must also be
present. The function reports an error rather than silently completing a
nonhierarchical parameter set.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34**, 44,
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9).
