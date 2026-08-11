# Find unique patterns in matrix columns

Given a matrix (for example of bootstrap replications) construct the
matrix of unique patterns of non-zeroes, together with a vector of
pointers back to that matrix.

## Usage

``` r
find_unique_patterns(x)
```

## Arguments

- x:

  a matrix

## Value

The original data `x` with the additional components

- yuniq:

  matrix of unique patterns of non-zeroes/zeroes in the columns of `x`

- pointers:

  vector of length dim(x)\[2\] giving the column of yuniq corresponding
  to each column of x

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
[\doi{10.1007/s11222-023-10346-9}](NA).
