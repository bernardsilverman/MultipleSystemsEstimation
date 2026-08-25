# Population estimation using BIC model selection

Selects a model using the BIC criterion, accounting for sparse data and
checking parameter identifiability and existence of the extended MLE.
Hierarchical log-linear models are fitted to the observed data and
ordered by increasing BIC. The point estimate is based on the model with
the lowest BIC value.

If requested, constructs bootstrap BCa confidence intervals for
population size that allow for uncertainty arising from BIC-based model
selection. To reduce the computational load to feasible levels, an
approximation is used in which the bootstrap and jackknife calculations
are restricted to a specified number of the best-ranked models.

## Usage

``` r
estimate_population_bic(
  zdat,
  nboot = 0,
  iseed = 1234,
  alpha = c(0.025, 0.1, 0.9, 0.975),
  maxorder = NULL,
  ntopmodels = NULL,
  return_details = FALSE
)
```

## Arguments

- zdat:

  A capture history data matrix with \\t+1\\ columns. The first \\t\\
  columns correspond to the capture lists and contain zeros and ones
  defining the observed capture histories. The final column contains the
  number of cases having each capture history. Capture histories not
  explicitly included in the data are assumed to have zero count.

- nboot:

  Non-negative integer giving the number of bootstrap replications. If
  `nboot = 0`, the best-BIC population estimate and model are returned
  together with the complete original-data BIC enumeration, but no
  bootstrap or jackknife calculations are performed.

- iseed:

  Integer seed used to initialise the random-number generator. The
  default is 1234.

- alpha:

  Numeric vector of cumulative probability levels at which the endpoints
  of the BCa confidence intervals are to be evaluated. The default is
  `c(0.025, 0.1, 0.9, 0.975)`.

- maxorder:

  Maximum order of interaction to include in the hierarchical models
  considered. If `NULL`, an automatic choice is made according to the
  number of lists: 1 for two lists, 2 for three lists, 3 for four or
  five lists, and 2 for six lists.

- ntopmodels:

  Number of models, ranked by BIC on the original data, to retain for
  bootstrap inference. If `NULL`, an automatic choice is made when
  `nboot > 0`: all available models for two or three lists, 20 models
  for four lists, and 100 models for five or six lists. This argument is
  irrelevant when `nboot = 0`, because no model subsetting is then
  performed. If `ntopmodels` is greater than the total number of
  available models, then all models are considered.

- return_details:

  Logical. If `TRUE`, include the bootstrap estimates, BCa acceleration,
  selected-model BIC, complete original-data BIC enumeration, and
  effects estimated at minus infinity. The default is `FALSE`.

## Value

A list with components:

- `input`:

  A list containing the original `call` and `data`.

- `method`:

  The character string `"bic"`.

- `estimate`:

  A named numeric vector containing the estimated `dark_figure` and
  `total` population.

- `fitted_model`:

  The hierarchy with the smallest BIC on the original data.

- `uncertainty`:

  A two-row matrix of BCa endpoints for the dark figure and total
  population when `nboot > 0`; otherwise an explanatory character
  string.

- `details`:

  If `return_details = TRUE`, a list containing
  `minus_infinity_effects`, `bootstrap_estimates`, `bca_acceleration`,
  `BIC`, and `bic_results`. The last is the complete original-data
  enumeration, with one row per eligible model giving its
  total-population estimate, BIC and maximum interaction order, ordered
  by increasing BIC. Otherwise `"not requested"`.

## Details

If `nboot > 0`, the routine implements the bootstrap procedure described
by Silverman, Chan and Vincent (2024). Multinomial bootstrap samples are
generated and the BIC model-selection procedure is repeated for each
bootstrap sample, restricting consideration to a set of models having
the smallest BIC values for the original data.

With the default settings, all available models are retained for
bootstrap inference for two- and three-list data. For four-list data,
the 20 models with the smallest BIC values are retained, and for five-
and six-list data the 100 models with the smallest BIC values are
retained. Focusing attention on a subset of all possible models enables
considerable computational economies.

The possible hierarchical models are drawn from an exhaustive model
catalogue within the package. This contains all hierarchical models for
up to five lists, but for six lists is restricted to models with
interactions of order at most 2. There are 32,768 such six-list models.
Allowing interactions of order 3 would increase this to 3,702,013
models, making exhaustive enumeration and fitting computationally
impractical.

For data with more than six lists, the routine stops with an informative
error.

A small value of `nboot`, such as that used in the example, is useful
only for checking that the routine runs. A substantially larger number
of bootstrap replications should be used for substantive inference.

## References

Silverman, B. W., Chan, L. and Vincent, K. (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection. *Statistics
and Computing*, **34**, 44. Available from
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9).

## Examples

``` r
data(Korea)

# A very small number of bootstrap replications is used here only
# to keep the example quick.
estimate_population_bic(Korea, nboot = 10)
#> $input
#> $input$call
#> estimate_population_bic(zdat = Korea, nboot = 10)
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
#>                 0.025       0.1      0.9 0.975
#> dark_figure  11.72727  11.82012 197.3872   213
#> total       134.72727 134.82012 320.3872   336
#> 
#> $details
#> [1] "not requested"
#> 
```
