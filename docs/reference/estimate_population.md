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

  Capture-pattern data. The first columns identify list membership and
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
  hiermod = "[12,23]"
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
#> $hiermod
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
