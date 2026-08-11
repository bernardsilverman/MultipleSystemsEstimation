# Search subsets for a property which is inherited in a particular way

The following routine returns a list of all subsets of a given set for
which a specified property is `TRUE`. It is assumed that if the property
is `FALSE` for any particular subset, it is also `FALSE` for all
supersets of that subset. This enables a branch search strategy to be
used to obviate the need to search over supersets of subsets already
eliminated from consideration. It is used within the hierarchical search
step of
[`checkallmodels`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/checkallmodels.md).

## Usage

``` r
subsetsearch(n, checkfun, testnull = TRUE, ...)
```

## Arguments

- n:

  an integer such that the search is over all subsets of \\\\1,...,n\\\\

- checkfun:

  a function which takes arguments `zset`, a subset of \\\\1,...,n\\\\
  and ... The function returns the value `TRUE` or `FALSE`. It needs to
  have the property that if it is `FALSE` for any particular `zset`, it
  is also `FALSE` for all supersets of `zset`.

- testnull:

  If `TRUE`, then include the null subset in the search. If `FALSE`, do
  not test the null subset.

- ...:

  other arguments

## Value

A list of all subsets for which the value is `TRUE`

## References

Chan, L., Silverman, B. W., and Vincent, K. (2021). Multiple Systems
Estimation for Sparse Capture Data: Inferential Challenges when there
are Non-Overlapping Lists. *Journal of the American Statistical
Association*, **116(535)**, 1297-1306, Available from
<https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748>.
