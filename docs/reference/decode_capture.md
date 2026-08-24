# Decode capture history

Converts an encoded capture history to a logical vector indicating
membership of each list.

## Usage

``` r
decode_capture(k, nlists)
```

## Arguments

- k:

  Integer encoding of a capture history.

- nlists:

  Number of lists.

## Value

A logical vector of length `nlists`.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34**, 44,
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9).
