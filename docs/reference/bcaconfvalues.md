# BCa confidence intervals

The BCa confidence intervals use percentiles of the bootstrap
distribution of the population size , but adjust the percentile actually
used. The adjusted percentiles depend on an estimated bias parameter,
and the quantile function of the estimated bias parameter is the
proportion of the bootstrap estimates that fall below the estimate from
the original data, and an estimated acceleration factor, which
derivation depends on a jackknife approach. This routine is called
internally by `estimatepopulation`.

## Usage

``` r
bcaconfvalues(
  bootreps,
  popest,
  ahat,
  alpha = c(0.025, 0.05, 0.1, 0.16, 0.84, 0.9, 0.95, 0.975)
)
```

## Arguments

- bootreps:

  Point estimates of total population sizes from each bootstrap sample.

- popest:

  A point estimate of the total population of the original data set.

- ahat:

  the estimated acceleration factor

- alpha:

  Bootstrap quantiles of interests

## Value

BCa confidence intervals

## References

Chan, L., Silverman, B. W., and Vincent, K. (2021). Multiple Systems
Estimation for Sparse Capture Data: Inferential Challenges when there
are Non-Overlapping Lists. *Journal of the American Statistical
Association*, **116(535)**, 1297-1306, Available from
<https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748>.

DiCiccio, T. J. and Efron, B. (1996). Bootstrap Confidence Intervals.
*Statistical Science*, **40(3)**, 189-228.

Efron, B. (1987). Better Bootstrap Confidence Intervals. *Journal of the
American Statistical Association*, **82(397)**, 171-185.
