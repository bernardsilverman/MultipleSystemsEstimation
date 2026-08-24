# Encode capture history

Encodes a binary capture history as the integer \$\$1 + \sum\_{i \in S}
2^{i-1},\$\$ where \\S\\ is the set of lists containing the case. Thus 1
represents the intercept or empty set, 2 and 3 represent lists 1 and 2
respectively, and 4 represents the two-list history 12.

## Usage

``` r
encode_capture(z)
```

## Arguments

- z:

  Logical vector, or vector of zeros and ones, defining a capture
  history.

## Value

The integer encoding of the capture history.
