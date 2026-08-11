# Decode capture history

Given a capture history as a number and the number of lists, decode it
into a logical vector giving presence or absence in the capture history.

## Usage

``` r
decode_capture(k, nlists)
```

## Arguments

- k:

  The capture history to be decoded

- nlists:

  The number of lists

## Value

A logical vector of length `nlists` giving presence or absence in the
capture history

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
[\doi{10.1007/s11222-023-10346-9}](NA).

## Examples

``` r
decode_capture(2,5)
#> [1]  TRUE FALSE FALSE FALSE FALSE
decode_capture(1,4)
#> [1] FALSE FALSE FALSE FALSE
```
