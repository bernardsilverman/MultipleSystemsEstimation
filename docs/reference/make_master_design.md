# Set up the inclusion matrix for all possible capture histories

This is the master design matrix which maps parameters to observations.
Rows correspond to observations and columns to parameters.

## Usage

``` r
make_master_design(nlists)
```

## Arguments

- nlists:

  Number of lists.

## Value

A binary matrix whose \\(i,j)\\ element is 1 when the expected log count
for history \\i\\ depends on parameter \\j\\; equivalently, when \\j\\
is an ancestor of \\i\\.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34**, 44,
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9).
