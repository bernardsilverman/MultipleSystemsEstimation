# Check a subset of the parameter set theta

This routine leaves out a particular set of parameters corresponding to
the two-list effects from the parameter set theta. For the resulting
model, it constructs the linear programming problem to check whether the
extended maximum likelihood estimates of the parameters exists. It is
called internally by `checkallmodels`.

## Usage

``` r
checkthetasubset(zset, amat, tvec, nlists)
```

## Arguments

- zset:

  set of indices that is not included, numbered among the two-list
  effects only

- amat:

  a design matrix

- tvec:

  vector of sufficient statistics

- nlists:

  number of lists in the original capture-recapture matrix

## Value

If the return result is `TRUE`, the linear program shows that the
extended maximum likelihood estimate does not exist. If the return
result is `FALSE`, the estimate exists.

## References

Chan, L., Silverman, B. W., and Vincent, K. (2021). Multiple Systems
Estimation for Sparse Capture Data: Inferential Challenges when there
are Non-Overlapping Lists. *Journal of the American Statistical
Association*, **116(535)**, 1297-1306, Available from
<https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748>.
