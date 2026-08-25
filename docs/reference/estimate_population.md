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

A list with components:

- `input`:

  A list containing the original `call` and `data`.

- `method`:

  The method actually used, including after `method = "auto"` has been
  resolved.

- `estimate`:

  A named numeric vector containing the estimated `dark_figure` and
  `total` population.

- `fitted_model`:

  The supplied, selected or retained model in hierarchy notation.

- `uncertainty`:

  Confidence or credible endpoints for the dark figure and total
  population, or an explanatory character string when bootstrap
  inference was not requested.

- `details`:

  Method-specific supporting output when requested, and
  `"not requested"` otherwise.

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
#> $input
#> $input$call
#> estimate_population(zdat = Korea)
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
#> [1] "bic"
#> 
#> $estimate
#> dark_figure       total 
#>    34.16667   157.16667 
#> 
#> $fitted_model
#> [1] "[12,23]"
#> 
#> $uncertainty
#> [1] "not calculated because nboot = 0"
#> 
#> $details
#> [1] "not requested"
#> 

# Pass BIC-specific arguments through ...
estimate_population(
  Korea,
  method = "bic",
  nboot = 100
)
#> $input
#> $input$call
#> estimate_population(zdat = Korea, method = "bic", nboot = 100)
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
#> [1] "bic"
#> 
#> $estimate
#> dark_figure       total 
#>    34.16667   157.16667 
#> 
#> $fitted_model
#> [1] "[12,23]"
#> 
#> $uncertainty
#>                 0.025       0.1      0.9    0.975
#> dark_figure   6.54714  13.09813 162.0565 294.0102
#> total       129.54714 136.09813 285.0565 417.0102
#> 
#> $details
#> [1] "not requested"
#> 

# Pass stepwise-specific arguments through ...
estimate_population(
  Korea,
  method = "stepwise",
  pthresh = 0.02
)
#> $input
#> $input$call
#> estimate_population(zdat = Korea, method = "stepwise", pthresh = 0.02)
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

# The bayesthresh method requires the suggested MCMCpack package.
if (requireNamespace("MCMCpack", quietly = TRUE)) {
  estimate_population(
    Western,
    method = "bayesthresh",
    burnin = 100,
    mcmc = 1000
  )
}
#> $input
#> $input$call
#> estimate_population(zdat = Western, method = "bayesthresh", burnin = 100, 
#>     mcmc = 1000)
#> 
#> $input$data
#>    A B C D E   n
#> 1  1 0 0 0 0  52
#> 2  0 1 0 0 0  90
#> 3  0 0 1 0 0 114
#> 4  0 0 0 1 0  45
#> 5  0 0 0 0 1  21
#> 6  1 0 1 0 0   4
#> 7  1 0 0 1 0   2
#> 8  1 0 0 0 1   5
#> 9  0 1 1 0 0   6
#> 10 0 1 0 1 0   1
#> 11 0 0 0 1 1   3
#> 12 1 0 1 0 1   1
#> 13 0 1 1 1 0   1
#> 
#> 
#> $method
#> [1] "bayesthresh"
#> 
#> $estimate
#> dark_figure       total 
#>    2207.763    2552.763 
#> 
#> $fitted_model
#> [1] "[15,45,2,3]"
#> 
#> $uncertainty
#>                0.025      0.1      0.9    0.975
#> dark_figure 1331.599 1510.708 2957.839 3384.523
#> total       1676.599 1855.708 3302.839 3729.523
#> 
#> $details
#> [1] "not requested"
#> 

# Pass a fixed-model specification through ...
estimate_population(
  Korea,
  method = "fixed",
  model = "[12,23]"
)
#> $input
#> $input$call
#> estimate_population(zdat = Korea, method = "fixed", model = "[12,23]")
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
#> [1] "fixed"
#> 
#> $estimate
#> dark_figure       total 
#>    34.16667   157.16667 
#> 
#> $fitted_model
#> [1] "[12,23]"
#> 
#> $uncertainty
#> [1] "not calculated because nboot = 0"
#> 
#> $details
#> [1] "not requested"
#> 
```
