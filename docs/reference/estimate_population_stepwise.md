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

  The model selected and fitted to the original data.

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

The stepwise procedure considers two-list interactions only;
higher-order interactions are not candidates for selection.

The procedure is first applied to the observed data to obtain the point
estimate and fitted model.

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
#> $MSEfit$fit$coefficients
#>         1         2         3         4         5 
#>  4.982083 -2.832024 -3.620482  5.412241 -1.268511 
#> 
#> $MSEfit$fit$residuals
#>             2             3             4             5             6 
#> -4.176136e-01  2.812500e-01  4.829545e-02 -1.733031e-16  1.484848e+00 
#>             7             8 
#> -1.000000e+00 -1.717172e-01 
#> 
#> $MSEfit$fit$fitted.values
#>         2         3         4         5         6         7         8 
#>  8.585366  3.902439 51.512195 41.000000  2.414634  1.097561 14.487805 
#> 
#> $MSEfit$fit$effects
#>           1           2           3           4           5             
#> -37.8739059  -0.1087053  -1.3877657  -7.6927803   4.7545965  -1.7833522 
#>             
#>  -2.3652222 
#> 
#> $MSEfit$fit$R
#>           1         2         3         4          5
#> 1 -11.09054 -6.942856 -6.401854 -5.951020 -5.3198508
#> 2   0.00000  5.366260  4.016362  4.599652 -3.7330505
#> 3   0.00000  0.000000 -3.726271 -2.530297  0.9334429
#> 4   0.00000  0.000000  0.000000 -1.739589 -1.3577277
#> 5   0.00000  0.000000  0.000000  0.000000 -3.7481703
#> 
#> $MSEfit$fit$rank
#> [1] 5
#> 
#> $MSEfit$fit$qr
#> $qr
#>              1           2          3           4          5
#> 2 -11.09053653 -6.94285621 -6.4018544 -5.95101958 -5.3198508
#> 3   0.17812116  5.36626016  4.0163624  4.59965154 -3.7330505
#> 4   0.64714630 -0.39565940 -3.7262707 -2.53029702  0.9334429
#> 5   0.57735027  0.84023239 -0.2576296 -1.73958866 -1.3577277
#> 6   0.14011129 -0.08566277 -0.2441369 -0.20535175 -3.7481703
#> 7   0.09446301  0.13747440  0.2389992 -0.29265174  0.3022582
#> 8   0.34320115 -0.20983008  0.4234623  0.06809711  0.9058907
#> 
#> $rank
#> [1] 5
#> 
#> $qraux
#> [1] 1.264196 1.259224 1.798488 1.931423 1.296651
#> 
#> $pivot
#> [1] 1 2 3 4 5
#> 
#> $tol
#> [1] 1e-11
#> 
#> attr(,"class")
#> [1] "qr"
#> 
#> $MSEfit$fit$family
#> 
#> Family: poisson 
#> Link function: log 
#> 
#> 
#> $MSEfit$fit$linear.predictors
#>          2          3          4          5          6          7          8 
#> 2.15005911 1.36160175 3.94181858 3.71357207 0.88154778 0.09309042 2.67330725 
#> 
#> $MSEfit$fit$deviance
#> [1] 8.566946
#> 
#> $MSEfit$fit$aic
#> [1] 44.90767
#> 
#> $MSEfit$fit$null.deviance
#> [1] 143.5474
#> 
#> $MSEfit$fit$iter
#> [1] 5
#> 
#> $MSEfit$fit$weights
#>         2         3         4         5         6         7         8 
#>  8.585366  3.902439 51.512195 41.000000  2.414634  1.097561 14.487805 
#> 
#> $MSEfit$fit$prior.weights
#> 2 3 4 5 6 7 8 
#> 1 1 1 1 1 1 1 
#> 
#> $MSEfit$fit$df.residual
#> [1] 2
#> 
#> $MSEfit$fit$df.null
#> [1] 6
#> 
#> $MSEfit$fit$y
#>  2  3  4  5  6  7  8 
#>  5  5 54 41  6  0 12 
#> 
#> $MSEfit$fit$converged
#> [1] TRUE
#> 
#> $MSEfit$fit$boundary
#> [1] FALSE
#> 
#> $MSEfit$fit$abundance
#>        1 
#> 268.7778 
#> 
#> $MSEfit$fit$bic
#> [1] 58.9686
#> 
#> $MSEfit$fit$neginfpars
#> numeric(0)
#> 
#> 
#> $MSEfit$hiermod
#> [1] "[12,3]"
#> 
#> $MSEfit$selected
#> [1]  TRUE FALSE FALSE
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
#> $MSEfit$fit$coefficients
#>         1         2         3         4         5 
#>  4.982083 -2.832024 -3.620482  5.412241 -1.268511 
#> 
#> $MSEfit$fit$residuals
#>             2             3             4             5             6 
#> -4.176136e-01  2.812500e-01  4.829545e-02 -1.733031e-16  1.484848e+00 
#>             7             8 
#> -1.000000e+00 -1.717172e-01 
#> 
#> $MSEfit$fit$fitted.values
#>         2         3         4         5         6         7         8 
#>  8.585366  3.902439 51.512195 41.000000  2.414634  1.097561 14.487805 
#> 
#> $MSEfit$fit$effects
#>           1           2           3           4           5             
#> -37.8739059  -0.1087053  -1.3877657  -7.6927803   4.7545965  -1.7833522 
#>             
#>  -2.3652222 
#> 
#> $MSEfit$fit$R
#>           1         2         3         4          5
#> 1 -11.09054 -6.942856 -6.401854 -5.951020 -5.3198508
#> 2   0.00000  5.366260  4.016362  4.599652 -3.7330505
#> 3   0.00000  0.000000 -3.726271 -2.530297  0.9334429
#> 4   0.00000  0.000000  0.000000 -1.739589 -1.3577277
#> 5   0.00000  0.000000  0.000000  0.000000 -3.7481703
#> 
#> $MSEfit$fit$rank
#> [1] 5
#> 
#> $MSEfit$fit$qr
#> $qr
#>              1           2          3           4          5
#> 2 -11.09053653 -6.94285621 -6.4018544 -5.95101958 -5.3198508
#> 3   0.17812116  5.36626016  4.0163624  4.59965154 -3.7330505
#> 4   0.64714630 -0.39565940 -3.7262707 -2.53029702  0.9334429
#> 5   0.57735027  0.84023239 -0.2576296 -1.73958866 -1.3577277
#> 6   0.14011129 -0.08566277 -0.2441369 -0.20535175 -3.7481703
#> 7   0.09446301  0.13747440  0.2389992 -0.29265174  0.3022582
#> 8   0.34320115 -0.20983008  0.4234623  0.06809711  0.9058907
#> 
#> $rank
#> [1] 5
#> 
#> $qraux
#> [1] 1.264196 1.259224 1.798488 1.931423 1.296651
#> 
#> $pivot
#> [1] 1 2 3 4 5
#> 
#> $tol
#> [1] 1e-11
#> 
#> attr(,"class")
#> [1] "qr"
#> 
#> $MSEfit$fit$family
#> 
#> Family: poisson 
#> Link function: log 
#> 
#> 
#> $MSEfit$fit$linear.predictors
#>          2          3          4          5          6          7          8 
#> 2.15005911 1.36160175 3.94181858 3.71357207 0.88154778 0.09309042 2.67330725 
#> 
#> $MSEfit$fit$deviance
#> [1] 8.566946
#> 
#> $MSEfit$fit$aic
#> [1] 44.90767
#> 
#> $MSEfit$fit$null.deviance
#> [1] 143.5474
#> 
#> $MSEfit$fit$iter
#> [1] 5
#> 
#> $MSEfit$fit$weights
#>         2         3         4         5         6         7         8 
#>  8.585366  3.902439 51.512195 41.000000  2.414634  1.097561 14.487805 
#> 
#> $MSEfit$fit$prior.weights
#> 2 3 4 5 6 7 8 
#> 1 1 1 1 1 1 1 
#> 
#> $MSEfit$fit$df.residual
#> [1] 2
#> 
#> $MSEfit$fit$df.null
#> [1] 6
#> 
#> $MSEfit$fit$y
#>  2  3  4  5  6  7  8 
#>  5  5 54 41  6  0 12 
#> 
#> $MSEfit$fit$converged
#> [1] TRUE
#> 
#> $MSEfit$fit$boundary
#> [1] FALSE
#> 
#> $MSEfit$fit$abundance
#>        1 
#> 268.7778 
#> 
#> $MSEfit$fit$bic
#> [1] 58.9686
#> 
#> $MSEfit$fit$neginfpars
#> numeric(0)
#> 
#> 
#> $MSEfit$hiermod
#> [1] "[12,3]"
#> 
#> $MSEfit$selected
#> [1]  TRUE FALSE FALSE
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
