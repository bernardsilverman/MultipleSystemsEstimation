# Estimate population size

Provides a common interface to the principal population-estimation
methods in MultipleSystemsEstimation.

## Usage

``` r
estimate_population(
  zdat,
  method = c("auto", "bic", "stepwise", "fixed", "bayesthresh"),
  ...
)
```

## Arguments

- zdat:

  Capture history data. The first columns identify list membership and
  the final column contains the observed counts.

- method:

  Estimation method. One of `"auto"`, `"bic"`, `"stepwise"`,
  `"bayesthresh"`, or `"fixed"`. The individual methods and their
  available arguments are documented in
  [`estimate_population_bic`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_bic.md),
  [`estimate_population_stepwise`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_stepwise.md),
  [`estimate_population_bayesthresh`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_bayesthresh.md),
  and
  [`estimate_population_fixed`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_fixed.md).

- ...:

  Additional arguments passed to the selected estimation function. See
  [`estimate_population_bic`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_bic.md),
  [`estimate_population_stepwise`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_stepwise.md),
  [`estimate_population_bayesthresh`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_bayesthresh.md),
  and
  [`estimate_population_fixed`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_fixed.md)
  for the arguments available for each method.

## Value

The object returned by
[`estimate_population_bic`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_bic.md),
[`estimate_population_stepwise`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_stepwise.md),
[`estimate_population_bayesthresh`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_bayesthresh.md),
or
[`estimate_population_fixed`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_fixed.md).

## Details

The default method is `"auto"`. This selects the BIC-based method when
there are no more than five lists, and the stepwise method when there
are six or more lists.

For six-list data, the BIC method remains available by specifying
`method = "bic"`, but exhaustive BIC enumeration is computationally
burdensome. Accordingly, `method = "auto"` selects the stepwise method
for six-list data and issues an informational message.

The estimation method can be selected explicitly using `method = "bic"`,
`"stepwise"`, `"bayesthresh"`, or `"fixed"`. The `"bayesthresh"` method
requires the suggested package MCMCpack.

Through their respective `maxorder` arguments, the stepwise, BIC and
Bayesian-threshold methods can all be restricted to two-list
interactions, or consider higher-order interactions.

Method-specific arguments are passed through `...` to the selected
estimation function. See the documentation for the individual methods
for details of the available arguments and their defaults.

## See also

[`estimate_population_bic`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_bic.md),
[`estimate_population_stepwise`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_stepwise.md),
[`estimate_population_bayesthresh`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_bayesthresh.md),
[`estimate_population_fixed`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_fixed.md)

## Examples

``` r
data(Korea)

# Three lists: automatically uses the BIC method.
estimate_population(Korea)
#> $popest
#> [1] 157.1667
#> 
#> $model
#> [1] "[12,23]"
#> 
#> $BIC
#> [1] 57.14081
#> 
#> $bic_results
#> $bic_results$res
#>         abundance       BIC modelsorder
#> [12,23]  157.1667  57.14081           2
#> [12,3]   268.7778  58.96860           2
#> [13,23]  123.4630  91.89711           2
#> [23,1]   126.1944 138.06368           2
#> [13,2]   126.9394 155.16018           2
#> [1,2,3]  141.9926 184.14431           1
#> 
#> $bic_results$xdata
#>      b c d Count
#> [1,] 0 0 1    41
#> [2,] 0 1 0     5
#> [3,] 1 0 0     5
#> [4,] 1 1 0    54
#> [5,] 1 0 1     6
#> [6,] 0 1 1     0
#> [7,] 1 1 1    12
#> 
#> $bic_results$maxorder
#> [1] 2
#> 
#> 
#> $BCaquantiles
#> NULL
#> 
#> attr(,"method")
#> [1] "bic"
#> attr(,"nlists")
#> [1] 3

# Pass BIC-specific arguments through ...
estimate_population(
  Korea,
  method = "bic",
  nboot = 100
)
#> $popest
#> [1] 157.1667
#> 
#> $model
#> [1] "[12,23]"
#> 
#> $BIC
#> [1] 57.14081
#> 
#> $bic_results
#> $bic_results$res
#>         abundance       BIC modelsorder
#> [12,23]  157.1667  57.14081           2
#> [12,3]   268.7778  58.96860           2
#> [13,23]  123.4630  91.89711           2
#> [23,1]   126.1944 138.06368           2
#> [13,2]   126.9394 155.16018           2
#> [1,2,3]  141.9926 184.14431           1
#> 
#> $bic_results$xdata
#>      b c d Count
#> [1,] 0 0 1    41
#> [2,] 0 1 0     5
#> [3,] 1 0 0     5
#> [4,] 1 1 0    54
#> [5,] 1 0 1     6
#> [6,] 0 1 1     0
#> [7,] 1 1 1    12
#> 
#> $bic_results$maxorder
#> [1] 2
#> 
#> 
#> $BCaquantiles
#>            0.025      0.1      0.9    0.975
#> [12,23] 133.2116 136.9897 199.7547 303.2850
#> [12,3]  133.1861 136.8837 282.6830 409.4789
#> [13,23] 133.1861 136.8837 282.6830 409.4789
#> [23,1]  129.5471 136.0981 285.0565 417.0102
#> [13,2]  129.5471 136.0981 285.0565 417.0102
#> [1,2,3] 129.5471 136.0981 285.0565 417.0102
#> 
#> attr(,"method")
#> [1] "bic"
#> attr(,"nlists")
#> [1] 3

# Pass stepwise-specific arguments through ...
estimate_population(
  Korea,
  method = "stepwise",
  pthresh = 0.02
)
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
#> attr(,"method")
#> [1] "stepwise"
#> attr(,"nlists")
#> [1] 3

# The bayesthresh method requires the suggested MCMCpack package.
if (requireNamespace("MCMCpack", quietly = TRUE)) {
  estimate_population(
    Western,
    method = "bayesthresh",
    burnin = 100,
    mcmc = 1000
  )
}
#> $call
#> estimate_population_bayesthresh(zdat = zdat, burnin = 100, mcmc = 1000)
#> 
#> $popest
#> [1] 2552.763
#> 
#> $quantiles
#>     2.5%      10%      50%      90%    97.5% 
#> 1676.599 1855.708 2552.763 3302.839 3729.523 
#> 
#> $retained_interactions
#> [1] "A:E" "D:E"
#> 
#> $threshold_statistics
#>       A:B       A:C       B:C       A:D       B:D       C:D       A:E       B:E 
#> 1.3057152 0.9920623 0.6130628 0.2974010 0.3391920 1.1302809 3.6826238 0.9032479 
#>       C:E       D:E 
#> 0.1200008 2.3012798 
#> 
#> attr(,"method")
#> [1] "bayesthresh"
#> attr(,"nlists")
#> [1] 5

# Pass a fixed-model specification through ...
estimate_population(
  Korea,
  method = "fixed",
  model = "[12,23]"
)
#> $popest
#> [1] 157.1667
#> 
#> $MSEfit
#> $MSEfit$coefficients
#>          1          2          3          4          5          7 
#>  3.5312505 -1.9218126 -2.1069550  4.5020294  0.1823216 -1.7749524 
#> 
#> $MSEfit$residuals
#>             2             3             4             5             6 
#>  1.243450e-15  2.033898e-01 -1.540832e-02 -1.039819e-15 -1.480297e-16 
#>             7             8 
#> -1.000000e+00  7.575758e-02 
#> 
#> $MSEfit$fitted.values
#>          2          3          4          5          6          7          8 
#>  5.0000000  4.1549296 54.8450704 41.0000000  6.0000000  0.8450704 11.1549296 
#> 
#> $MSEfit$effects
#>            1            2            3            4            5            7 
#> -38.17324951   0.05565632  -2.06756935  -7.97592824   4.31707491   2.59747923 
#>              
#>   1.04592269 
#> 
#> $MSEfit$R
#>           1         2         3         4          5             7
#> 1 -11.09054 -6.942852 -6.401866 -5.951016 -5.3198514 -1.082007e+00
#> 2   0.00000  5.366267  4.016349  4.599657 -3.6859875  6.788207e-01
#> 3   0.00000  0.000000 -3.726289 -2.530285  1.9463634 -6.297963e-01
#> 4   0.00000  0.000000  0.000000 -1.739607 -0.7907304 -2.931683e-16
#> 5   0.00000  0.000000  0.000000  0.000000 -3.5635856 -2.798269e+00
#> 7   0.00000  0.000000  0.000000  0.000000  0.0000000 -1.463408e+00
#> 
#> $MSEfit$rank
#> [1] 6
#> 
#> $MSEfit$qr
#> $qr
#>              1          2          3           4          5             7
#> 2 -11.09054482 -6.9428521 -6.4018662 -5.95101620 -5.3198514 -1.082007e+00
#> 3   0.18379613  5.3662666  4.0163494  4.59965721 -3.6859875  6.788207e-01
#> 4   0.66775306 -0.4295199 -3.7262891 -2.53028453  1.9463634 -6.297963e-01
#> 5   0.57734984  0.8218480 -0.3085109 -1.73960681 -0.7907304 -2.931683e-16
#> 6   0.22086289 -0.1420660 -0.4051436 -0.35737012 -3.5635856 -2.798269e+00
#> 7   0.08288986  0.1179924  0.2024121 -0.25354332  0.2840594 -1.463408e+00
#> 8   0.30114865 -0.1937084  0.3438908  0.04316316  0.8884669  9.534136e-01
#> 
#> $rank
#> [1] 6
#> 
#> $qraux
#> [1] 1.201619 1.261631 1.762527 1.897853 1.360468 1.301666
#> 
#> $pivot
#> [1] 1 2 3 4 5 6
#> 
#> $tol
#> [1] 1e-11
#> 
#> attr(,"class")
#> [1] "qr"
#> 
#> $MSEfit$family
#> 
#> Family: poisson 
#> Link function: log 
#> 
#> 
#> $MSEfit$linear.predictors
#>          2          3          4          5          6          7          8 
#>  1.6094379  1.4242955  4.0045123  3.7135721  1.7917595 -0.1683353  2.4118815 
#> 
#> $MSEfit$deviance
#> [1] 1.926975
#> 
#> $MSEfit$aic
#> [1] 40.2677
#> 
#> $MSEfit$null.deviance
#> [1] 143.5474
#> 
#> $MSEfit$iter
#> [1] 4
#> 
#> $MSEfit$weights
#>          2          3          4          5          6          7          8 
#>  5.0000000  4.1550716 54.8450625 41.0000000  6.0000000  0.8451009 11.1549494 
#> 
#> $MSEfit$prior.weights
#> 2 3 4 5 6 7 8 
#> 1 1 1 1 1 1 1 
#> 
#> $MSEfit$df.residual
#> [1] 1
#> 
#> $MSEfit$df.null
#> [1] 6
#> 
#> $MSEfit$y
#>  2  3  4  5  6  7  8 
#>  5  5 54 41  6  0 12 
#> 
#> $MSEfit$converged
#> [1] TRUE
#> 
#> $MSEfit$boundary
#> [1] FALSE
#> 
#> $MSEfit$abundance
#>        1 
#> 157.1667 
#> 
#> $MSEfit$bic
#> [1] 57.14081
#> 
#> $MSEfit$neginfpars
#> numeric(0)
#> 
#> 
#> $model
#> [1] "[12,23]"
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
#> attr(,"method")
#> [1] "fixed"
#> attr(,"nlists")
#> [1] 3
```
