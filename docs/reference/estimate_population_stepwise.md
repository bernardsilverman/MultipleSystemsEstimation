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
  alpha = c(0.025, 0.1, 0.9, 0.975),
  return_details = FALSE
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

- return_details:

  Logical. If `TRUE`, include the bootstrap estimates, BCa acceleration,
  and effects estimated at minus infinity. The default is `FALSE`.

## Value

A list with components:

- `input`:

  A list containing the original `call` and `data`.

- `method`:

  The character string `"stepwise"`.

- `estimate`:

  A named numeric vector containing the estimated `dark_figure` and
  `total` population.

- `fitted_model`:

  The hierarchy selected from the original data.

- `uncertainty`:

  A two-row matrix of BCa endpoints for the dark figure and total
  population when `nboot > 0`; otherwise an explanatory character
  string.

- `details`:

  If `return_details = TRUE`, a list containing
  `minus_infinity_effects`, `bootstrap_estimates`, and
  `bca_acceleration`. Otherwise `"not requested"`.

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
#> $input
#> $input$call
#> estimate_population_stepwise(zdat = Korea)
#> 
#> $input$data
#>      b c d Count
#> [1,] 0 0 1    41
#> [2,] 0 1 0     5
#> [3,] 1 0 0     5
#> [4,] 1 1 0    54
#> [5,] 1 0 1     6
#> [6,] 0 1 1     0
#> [7,] 1 1 1    12
#> 
#> 
#> $method
#> [1] "stepwise"
#> 
#> $estimate
#> dark_figure       total 
#>    145.7778    268.7778 
#> 
#> $fitted_model
#> [1] "[12,3]"
#> 
#> $uncertainty
#> [1] "not calculated because nboot = 0"
#> 
#> $details
#> [1] "not requested"
#> 

# Allow interactions of any order
estimate_population_stepwise(Kosovo, maxorder = Inf)
#> $input
#> $input$call
#> estimate_population_stepwise(zdat = Kosovo, maxorder = Inf)
#> 
#> $input$data
#>    EXH ABA OSCE HRW Frequency
#> 1    1   0    0   0      1131
#> 2    0   1    0   0       845
#> 3    0   0    1   0       936
#> 4    0   0    0   1       306
#> 5    1   1    0   0       177
#> 6    1   0    1   0       228
#> 7    1   0    0   1       106
#> 8    0   1    1   0       217
#> 9    0   1    0   1        31
#> 10   0   0    1   1       123
#> 11   1   1    1   0       181
#> 12   1   1    0   1        18
#> 13   1   0    1   1        42
#> 14   0   1    1   1        32
#> 15   1   1    1   1        27
#> 
#> 
#> $method
#> [1] "stepwise"
#> 
#> $estimate
#> dark_figure       total 
#>    13993.31    18393.31 
#> 
#> $fitted_model
#> [1] "[134,12,23]"
#> 
#> $uncertainty
#> [1] "not calculated because nboot = 0"
#> 
#> $details
#> [1] "not requested"
#> 

# A very small number of bootstrap replications is used here only
# to keep the example quick.
estimate_population_stepwise(Korea, nboot = 10)
#> $input
#> $input$call
#> estimate_population_stepwise(zdat = Korea, nboot = 10)
#> 
#> $input$data
#>      b c d Count
#> [1,] 0 0 1    41
#> [2,] 0 1 0     5
#> [3,] 1 0 0     5
#> [4,] 1 1 0    54
#> [5,] 1 0 1     6
#> [6,] 0 1 1     0
#> [7,] 1 1 1    12
#> 
#> 
#> $method
#> [1] "stepwise"
#> 
#> $estimate
#> dark_figure       total 
#>    145.7778    268.7778 
#> 
#> $fitted_model
#> [1] "[12,3]"
#> 
#> $uncertainty
#>                  0.025      0.1      0.9    0.975
#> dark_figure   4.982845 113.3239 202.7143 202.7143
#> total       127.982845 236.3239 325.7143 325.7143
#> 
#> $details
#> [1] "not requested"
#> 
```
