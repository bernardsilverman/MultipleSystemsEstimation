# Population estimation using stepwise model selection

Estimates the total population, including the unobserved population,
using an extension of the stepwise model-selection procedure of Chan,
Silverman and Vincent (2021). Optional bootstrap and jackknife
calculations provide BCa confidence intervals while repeating the
stepwise selection procedure for each resampled data set.

## Usage

``` r
estimate_population_stepwise(
  zdat,
  nboot = 0,
  pthresh = 0.02,
  maxorder = 2,
  iseed = 1234,
  alpha = c(0.025, 0.1, 0.9, 0.975)
)
```

## Arguments

- zdat:

  A capture history data matrix with \\t+1\\ columns. The first \\t\\
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

- maxorder:

  Maximum order of interactions considered for selection. An integer of
  at least 2, or `Inf` to allow interactions of any order. The default
  is 2.

- iseed:

  Integer seed used to initialise the random-number generator when
  `nboot > 0`. The default is 1234.

- alpha:

  Numeric vector of cumulative probability levels at which the endpoints
  of the BCa confidence intervals are to be calculated. This argument is
  used only when `nboot > 0`. The default is
  `c(0.025, 0.1, 0.9, 0.975)`.

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

  The endpoints of the BCa confidence intervals at the cumulative
  probability levels specified by `alpha`. This is `NULL` when
  `nboot = 0`.

## Details

By default, the procedure considers two-list interactions only,
reproducing the method of Chan, Silverman and Vincent (2021).
Higher-order interactions may be considered by increasing `maxorder`, or
by setting `maxorder = Inf`. An interaction is considered only when all
its lower-order terms are already present, so every candidate model is
hierarchical.

The procedure starts with a main-effects model. At each stage it
considers eligible interactions not already included and adds the
interaction having the smallest p-value if that p-value does not exceed
`pthresh`. It stops when no eligible interaction meets the threshold.
Candidate models are checked for parameter identifiability and existence
of the extended MLE.

The procedure is first applied to the observed data to obtain the point
estimate and fitted model.

If `nboot > 0`, multinomial bootstrap samples are generated from the
observed capture history counts. The complete stepwise model-selection
procedure is repeated for each bootstrap sample, so the resulting
inference allows for variation in the selected model rather than
treating the model selected from the original data as fixed.

A delete-one jackknife calculation is also carried out to estimate the
acceleration parameter required for the BCa confidence intervals. The
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

# Allow interactions of any order
estimate_population_stepwise(Kosovo, maxorder = Inf)
#> $popest
#> [1] 18393.31
#> 
#> $MSEfit
#> $MSEfit$fit
#> $MSEfit$fit$coefficients
#>         1         2         3         4         5         6         7         9 
#>  9.546334 -2.526055 -2.806998  1.027940 -2.718600  1.183653  1.416918 -3.784884 
#>        10        13        14 
#>  1.428911  1.778187 -1.201823 
#> 
#> $MSEfit$fit$residuals
#>             2             3             4             5             6 
#>  1.063424e-02 -3.632597e-15 -6.300011e-02  1.397745e-02 -5.444914e-02 
#>             7             8             9            10            11 
#> -5.612183e-02  7.821067e-02 -3.715691e-02 -8.669459e-04  6.153231e-01 
#>            12            13            14            15            16 
#>  5.136022e-03 -8.813974e-03  3.246076e-02  3.538959e-02 -4.662660e-02 
#> 
#> $MSEfit$fit$fitted.values
#>          2          3          4          5          6          7          8 
#> 1119.09923  845.00000  188.90077  923.09745  241.12928  229.90255  167.87072 
#>          9         10         11         12         13         14         15 
#>  317.80879  106.09198   19.19121   17.90802  124.09376   40.67951   30.90624 
#>         16 
#>   28.32049 
#> 
#> $MSEfit$fit$effects
#>           1           2           3           4           5           6 
#> -417.000667   -5.661547   17.877238   21.526333   16.915394  -13.118303 
#>           7           9          10          13          14             
#>   11.337271   45.411463   -3.721823   10.732979   -5.382643    2.002626 
#>                                     
#>    1.981420    1.363165    1.259394 
#> 
#> $MSEfit$fit$R
#>           1         2          3          4          5          6            7
#> 1  -66.3325 -28.79433 -23.035467  -6.075454 -26.924963  -7.206121  -6.88953423
#> 2    0.0000  32.87684  -7.917152   6.936840  -9.042426   8.227816  -0.06656175
#> 3    0.0000   0.00000 -30.572634 -10.400449   7.680704  -3.118335  -9.73973587
#> 4    0.0000   0.00000   0.000000 -14.484467 -12.096992  -4.342834  -3.69349232
#> 5    0.0000   0.00000   0.000000   0.000000 -27.819965 -11.854460 -10.82048606
#> 6    0.0000   0.00000   0.000000   0.000000   0.000000  13.757269  -2.00555417
#> 7    0.0000   0.00000   0.000000   0.000000   0.000000   0.000000 -13.41344035
#> 9    0.0000   0.00000   0.000000   0.000000   0.000000   0.000000   0.00000000
#> 10   0.0000   0.00000   0.000000   0.000000   0.000000   0.000000   0.00000000
#> 13   0.0000   0.00000   0.000000   0.000000   0.000000   0.000000   0.00000000
#> 14   0.0000   0.00000   0.000000   0.000000   0.000000   0.000000   0.00000000
#>             9            10            13            14
#> 1  -10.326764 -2.909584e+00 -3.376927e+00 -1.040214e+00
#> 2   -3.174036  3.322110e+00 -8.588528e-01  1.187697e+00
#> 3    5.452104 -1.801128e-01  8.295664e-01 -4.501362e-01
#> 4   -4.295005 -2.508390e-01 -1.545774e+00 -6.268944e-01
#> 5    6.347289 -6.847051e-01 -3.603145e+00 -1.711209e+00
#> 6    6.853987  7.946100e-01  3.556245e-01  1.985882e+00
#> 7   -8.016895  2.220446e-16 -2.442491e-15 -6.661338e-16
#> 9  -19.198310 -8.974869e+00 -1.019215e+01 -3.075258e+00
#> 10   0.000000 -9.578833e+00  3.385949e+00 -3.282207e+00
#> 13   0.000000  0.000000e+00  8.962560e+00  4.330197e+00
#> 14   0.000000  0.000000e+00  0.000000e+00  4.478733e+00
#> 
#> $MSEfit$fit$rank
#> [1] 11
#> 
#> $MSEfit$fit$qr
#> $qr
#>               1            2             3            4             5
#> 2  -66.33249584 -28.79433342 -23.035466767  -6.07545361 -2.692496e+01
#> 3    0.43822991  32.87683630  -7.917151969   6.93684035 -9.042426e+00
#> 4    0.20720037  -0.15726522 -30.572634481 -10.40044851  7.680704e+00
#> 5    0.45803370   0.57648357  -0.597354476 -14.48446691 -1.209699e+01
#> 6    0.23409850  -0.17768092  -0.003847799  -0.12999099 -2.781997e+01
#> 7    0.22858386   0.28769682   0.197838761  -0.16407887  5.445697e-01
#> 8    0.19532649  -0.14825294   0.420583295   0.33234063  2.610608e-01
#> 9    0.26875504   0.33825647  -0.350502647   0.43135432 -5.412339e-01
#> 10   0.15527979  -0.11785747  -0.002552282  -0.08622427 -4.521807e-02
#> 11   0.06604270   0.08312168   0.057159791  -0.04740585 -1.311509e-04
#> 12   0.06379659  -0.04842166   0.137368887   0.10854749 -6.684691e-02
#> 13   0.16793787   0.21136746  -0.219019769   0.26954182  6.221975e-02
#> 14   0.09615271  -0.07297997  -0.001580430  -0.05339199  2.012615e-01
#> 15   0.08381018   0.10548392   0.072537499  -0.06015945  1.996662e-01
#> 16   0.08022764  -0.06089284   0.172748749   0.13650430  1.072271e-01
#>                6            7           9            10            13
#> 2   -7.206121144  -6.88953423 -10.3267636 -2.909584e+00 -3.376927e+00
#> 3    8.227815581  -0.06656175  -3.1740363  3.322110e+00 -8.588528e-01
#> 4   -3.118334743  -9.73973587   5.4521041 -1.801128e-01  8.295664e-01
#> 5   -4.342833516  -3.69349232  -4.2950045 -2.508390e-01 -1.545774e+00
#> 6  -11.854459554 -10.82048606   6.3472885 -6.847051e-01 -3.603145e+00
#> 7   13.757268763  -2.00555417   6.8539875  7.946100e-01  3.556245e-01
#> 8   -0.364132391 -13.41344035  -8.0168953  2.220446e-16 -2.442491e-15
#> 9   -0.647835190   0.62862258 -19.1983103 -8.974869e+00 -1.019215e+01
#> 10   0.031328235  -0.03652234   0.4790859 -9.578833e+00  3.385949e+00
#> 11  -0.027055013  -0.01258160   0.2151278 -1.407269e-01  8.962560e+00
#> 12  -0.001027648  -0.08052081   0.2415994  2.396815e-01  8.014457e-02
#> 13   0.094551055  -0.09389834   0.5382917 -3.830500e-01 -5.324370e-01
#> 14  -0.158302235  -0.03521137   0.3374943  3.847256e-01 -5.208399e-01
#> 15   0.214877312   0.15560002   0.3188182 -1.969222e-01 -2.352982e-01
#> 16  -0.149562324   0.28497484   0.2237992  3.765565e-01 -5.603308e-01
#>               14
#> 2  -1.040214e+00
#> 3   1.187697e+00
#> 4  -4.501362e-01
#> 5  -6.268944e-01
#> 6  -1.711209e+00
#> 7   1.985882e+00
#> 8  -6.661338e-16
#> 9  -3.075258e+00
#> 10 -3.282207e+00
#> 11  4.330197e+00
#> 12  4.478733e+00
#> 13  6.148003e-01
#> 14 -5.904357e-01
#> 15  2.343869e-01
#> 16 -4.388493e-01
#> 
#> $rank
#> [1] 11
#> 
#> $qraux
#>  [1] 1.504322 1.551558 1.446151 1.735148 1.490002 1.586056 1.693815 1.332318
#>  [9] 1.668905 1.263586 1.160874
#> 
#> $pivot
#>  [1]  1  2  3  4  5  6  7  8  9 10 11
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
#>        2        3        4        5        6        7        8        9 
#> 7.020279 6.739337 5.241222 6.827735 5.485333 5.437656 5.123194 5.761450 
#>       10       11       12       13       14       15       16 
#> 4.664306 2.954452 2.885249 4.821037 3.705725 3.430958 3.343586 
#> 
#> $MSEfit$fit$deviance
#> [1] 10.25032
#> 
#> $MSEfit$fit$aic
#> [1] 133.6743
#> 
#> $MSEfit$fit$null.deviance
#> [1] 5336.116
#> 
#> $MSEfit$fit$iter
#> [1] 4
#> 
#> $MSEfit$fit$weights
#>          2          3          4          5          6          7          8 
#> 1119.09923  845.00000  188.90077  923.09745  241.12928  229.90255  167.87072 
#>          9         10         11         12         13         14         15 
#>  317.80879  106.09198   19.19121   17.90802  124.09376   40.67951   30.90624 
#>         16 
#>   28.32049 
#> 
#> $MSEfit$fit$prior.weights
#>  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 
#>  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1 
#> 
#> $MSEfit$fit$df.residual
#> [1] 4
#> 
#> $MSEfit$fit$df.null
#> [1] 14
#> 
#> $MSEfit$fit$y
#>    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
#> 1131  845  177  936  228  217  181  306  106   31   18  123   42   32   27 
#> 
#> $MSEfit$fit$converged
#> [1] TRUE
#> 
#> $MSEfit$fit$boundary
#> [1] FALSE
#> 
#> $MSEfit$fit$abundance
#>        1 
#> 18393.31 
#> 
#> $MSEfit$fit$bic
#> [1] 203.9573
#> 
#> $MSEfit$fit$neginfpars
#> numeric(0)
#> 
#> 
#> $MSEfit$hiermod
#> [1] "[134,12,23]"
#> 
#> $MSEfit$selected
#>  [1]  TRUE  TRUE  TRUE  TRUE FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE
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
