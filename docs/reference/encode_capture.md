# Encode capture history

Given a 0/1 capture history \\S\\, encode it as \$\$1 + \sum\_{i \in S}
2^{i-1},\$\$ where \\S\\ is the set of list numbers. Thus 1 represents
the intercept (the empty set), 2 and 3 represent the single list capture
histories 1 and 2, and 4 represents the two-list history 12.

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
#> Error in encode_capture(c(1, 0, 0, 0, 0)): could not find function "encode_capture"
encode_capture(c(1,1,1,1,0))
#> Error in encode_capture(c(1, 1, 1, 1, 0)): could not find function "encode_capture"
encode_capture(c(TRUE,FALSE,TRUE,FALSE))
#> Error in encode_capture(c(TRUE, FALSE, TRUE, FALSE)): could not find function "encode_capture"
```
