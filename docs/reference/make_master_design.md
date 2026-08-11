# Set up the inclusion matrix for all possible capture histories

This is the master design matrix which maps parameters to observations.
Rows correspond to observations and columns to parameters.

## Usage

``` r
make_master_design(nlists)
```

## Arguments

- nlists:

  The number of lists

## Value

A matrix whose \\(i,j)\\ element is 1 if the expected log of observation
\\i\\ depends on parameter \\j\\, in other words if \\j\\ is an ancestor
of \\i\\.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
<https://doi.org/10.1007/s11222-023-10346-9>.
