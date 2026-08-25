# Find unique support patterns

Finds the distinct zero and nonzero patterns among the columns of a
matrix, together with indices mapping the original columns to those
patterns.

## Usage

``` r
find_unique_patterns(x)
```

## Arguments

- x:

  Numeric matrix whose columns are data vectors, typically bootstrap
  replications.

## Value

A list with components:

- `x`:

  The original matrix.

- `yuniq`:

  A binary matrix containing the distinct support patterns as columns.

- `pointers`:

  An integer vector mapping each column of `x` to the corresponding
  column of `yuniq`.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34**, 44,
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9)
.
