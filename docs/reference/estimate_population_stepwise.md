# Population estimation using stepwise model selection

Estimates the total population, including the unobserved population,
using the stepwise model-selection procedure of Chan, Silverman and
Vincent (2021). Optional bootstrap and jackknife calculations provide
BCa confidence limits while repeating the stepwise selection procedure
for each resampled data set.

## Usage

``` r
estimate_population_stepwise(
  zdat,
  nboot = 0,
  pthresh = 0.02,
  iseed = 1234,
  alpha = c(0.025, 0.1, 0.9, 0.975)
)
```

## Arguments

- zdat:

  A capture-history data matrix with \\t+1\\ columns. The first \\t\\
  columns correspond to the capture lists and contain zeros and ones
  defining the observed capture histories. The final column contains the
  number of cases having each capture history. List names `A`, `B`, and
  so on are constructed if they are not supplied. Capture histories not
  explicitly included in the data are assumed to have zero count.

- nboot:

  Non-negative integer giving the number of bootstrap replications. If
  `nboot = 0`, only the point estimate and fitted model are returned and
  no bootstrap or jackknife calculations are performed. The default is
  0.

- pthresh:

  P-value threshold used by the stepwise model-selection procedure. The
  default is 0.02.

- iseed:

  Integer seed used to initialise the random-number generator when
  `nboot > 0`. The default is 1234.

- alpha:

  Numeric vector of cumulative probability levels at which BCa
  confidence limits are to be calculated. This argument is used only
  when `nboot > 0`. The default is `c(0.025, 0.1, 0.9, 0.975)`.

## Value

A list with the following components:

- `popest`:

  The estimated total population for the original data, including the
  estimated unobserved population.

- `MSEfit`:

  The model selected and fitted to the original data, in the format
  returned by
  [`stepwisefit`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/stepwisefit.md).

- `bootreps`:

  A numeric vector containing the estimated total population from each
  bootstrap sample. This is `NULL` when `nboot = 0`.

- `ahat`:

  The estimated BCa acceleration parameter. This is `NULL` when
  `nboot = 0`.

- `BCaquantiles`:

  The BCa confidence limits at the cumulative probability levels
  specified by `alpha`. This is `NULL` when `nboot = 0`.

## Details

The stepwise model-selection procedure is first applied to the observed
data to obtain the point estimate and fitted model.

If `nboot > 0`, multinomial bootstrap samples are generated from the
observed capture-history counts. The complete stepwise model-selection
procedure is repeated for each bootstrap sample, so the resulting
inference allows for variation in the selected model rather than
treating the model selected from the original data as fixed.

A delete-one jackknife calculation is also carried out to estimate the
acceleration parameter required for the BCa confidence limits. The
jackknife calculation takes account of the number of individuals having
each observed capture history.

A small positive value of `nboot`, such as that used in the example, is
useful only for checking that the routine runs. A substantially larger
number of bootstrap replications should be used for substantive
inference.

## References

Chan, L., Silverman, B. W., and Vincent, K. (2021). Multiple Systems
Estimation for Sparse Capture Data: Inferential Challenges when there
are Non-Overlapping Lists. *Journal of the American Statistical
Association*, **116**(535), 1297–1306. Available from
<https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748>.

DiCiccio, T. J. and Efron, B. (1996). Bootstrap Confidence Intervals.
*Statistical Science*, **11**(3), 189–228.

## Examples

``` r
data(Korea)

# Point estimate and fitted model without bootstrapping
estimate_population_stepwise(Korea)
#> $popest
#> [1] 268.7778
#> 
#> $MSEfit
#> $MSEfit$fit
#> 
#> Call:  glm(formula = zz$modelform, family = poisson, data = zz$datamatrix, 
#>     x = TRUE)
#> 
#> Coefficients:
#> (Intercept)            b            c            d          b:c  
#>       4.982       -2.832       -3.620       -1.269        5.412  
#> 
#> Degrees of Freedom: 6 Total (i.e. Null);  2 Residual
#> Null Deviance:       143.5 
#> Residual Deviance: 8.567     AIC: 44.91
#> 
#> $MSEfit$emptyoverlaps
#> 
#> 
#> 
#> 
#> $MSEfit$poisspempty
#> NULL
#> 
#> 
#> $bootreps
#> NULL
#> 
#> $ahat
#> NULL
#> 
#> $BCaquantiles
#> NULL
#> 

# A very small number of bootstrap replications is used here only
# to keep the example quick.
estimate_population_stepwise(Korea, nboot = 10)
#> $popest
#> [1] 268.7778
#> 
#> $MSEfit
#> $MSEfit$fit
#> 
#> Call:  glm(formula = zz$modelform, family = poisson, data = zz$datamatrix, 
#>     x = TRUE)
#> 
#> Coefficients:
#> (Intercept)            b            c            d          b:c  
#>       4.982       -2.832       -3.620       -1.269        5.412  
#> 
#> Degrees of Freedom: 6 Total (i.e. Null);  2 Residual
#> Null Deviance:       143.5 
#> Residual Deviance: 8.567     AIC: 44.91
#> 
#> $MSEfit$emptyoverlaps
#> 
#> 
#> 
#> 
#> $MSEfit$poisspempty
#> NULL
#> 
#> 
#> $bootreps
#>  [1] 231.8182 247.8000 288.7500 325.7143 240.0000 123.3125 295.0000 251.1000
#>  [9] 232.9091 123.5263
#> 
#> $ahat
#> [1] -0.02877338
#> 
#> $BCaquantiles
#>    0.025      0.1      0.9    0.975 
#> 127.9828 236.3239 325.7143 325.7143 
#> 
```
