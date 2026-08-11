# Encode capture history

Given a 0/1 capture history, encode it as number that corresponds to the
row number of the capture history data set

## Usage

``` r
encode_capture(z)
```

## Arguments

- z:

  The capture history to be encoded, as a logical vector or a vector of
  0s and 1s

## Value

The capture history encoded as a number that corresponds to the row
number of the capture history data set

## Examples

``` r
encode_capture(c(1,0,0,0,0))
#> [1] 2
encode_capture(c(1,1,1,1,0))
#> [1] 16
encode_capture(c(TRUE,FALSE,TRUE,FALSE))
#> [1] 6
```
