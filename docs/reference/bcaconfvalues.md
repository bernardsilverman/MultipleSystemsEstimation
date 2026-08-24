# BCa confidence intervals

Calculates bias-corrected and accelerated (BCa) confidence intervals
from bootstrap estimates of population size. The percentile levels are
adjusted using the proportion of bootstrap estimates below the
original-data estimate and an acceleration parameter obtained by
jackknife.

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

  Numeric vector of population estimates from the bootstrap samples.

- popest:

  Population estimate from the original data.

- ahat:

  Estimated BCa acceleration parameter.

- alpha:

  Cumulative probability levels at which the interval endpoints are
  required.

## Value

A named numeric vector containing the requested BCa confidence interval
endpoints.

## References

Chan, L., Silverman, B. W., and Vincent, K. (2021). Multiple Systems
Estimation for Sparse Capture Data: Inferential Challenges When There
Are Nonoverlapping Lists. *Journal of the American Statistical
Association*, **116**(535), 1297–1306.
[doi:10.1080/01621459.2019.1708748](https://doi.org/10.1080/01621459.2019.1708748).

DiCiccio, T. J. and Efron, B. (1996). Bootstrap confidence intervals.
*Statistical Science*, **11**(3), 189–228.

Efron, B. (1987). Better bootstrap confidence intervals. *Journal of the
American Statistical Association*, **82**(397), 171–185.
