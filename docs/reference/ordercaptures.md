# Order capture histories

Given a matrix with capture histories only, the routine orders the
capture histories first by the number of 1s in the capture history and
then lexicographically by columns.

## Usage

``` r
ordercaptures(zmat)
```

## Arguments

- zmat:

  Data matrix with \\t\\ columns. The \\t\\ columns, each corresponding
  to a particular list, are 0s and 1s defining the capture histories
  observed. Where a capture history is not explicitly listed, it is
  assumed that it has observed count zero.

## Value

A data matrix that is ordered first by the number of 1s in the capture
history and then lexicographically by columns.

## References

Chan, L., Silverman, B. W., and Vincent, K. (2021). Multiple Systems
Estimation for Sparse Capture Data: Inferential Challenges when there
are Non-Overlapping Lists. *Journal of the American Statistical
Association*, **116(535)**, 1297-1306, Available from
<https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748>.
