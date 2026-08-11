# Population estimation using BIC model selection

Constructs BCa confidence limits for population size while allowing for
uncertainty arising from BIC-based model selection.

## Usage

``` r
estimate_population_bic(
  zdat,
  nboot = 0,
  iseed = 1234,
  alpha = c(0.025, 0.1, 0.9, 0.975),
  maxorder = NULL,
  ntopmodels = NULL
)
```

## Arguments

- zdat:

  A capture-history data matrix with \\t+1\\ columns. The first \\t\\
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

  Numeric vector of cumulative probability levels at which the BCa
  confidence limits are to be evaluated. The default is
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
  performed.

## Value

A list with components:

- `popest`:

  The estimated total population for the original data, including the
  estimated unobserved population, from the model with the smallest BIC.

- `model`:

  The hierarchical model with the smallest BIC on the original data.

- `BIC`:

  The BIC value of the selected model.

- `bic_results`:

  The complete result of the original-data BIC enumeration, including
  all models allowed by `maxorder`, ordered by BIC.

- `BCaquantiles`:

  A matrix of BCa confidence limits for nested sets of the retained
  top-BIC models when `nboot > 0`, and `NULL` when `nboot = 0`.

When `nboot > 0`, the `BCaquantiles` component gives BCa inference based
on repeated BIC model selection. This is a numeric matrix of BCa
confidence limits. The columns correspond to the cumulative probability
levels supplied in `alpha`. The retained models are ordered by
increasing BIC for the original data. Row \\k\\ gives the confidence
limits obtained when model selection within each bootstrap replication
is restricted to the first \\k\\ models in this ordering. Thus, the
first row uses only the best-BIC model, the second row allows selection
between the best two models, and the final row allows selection among
all retained models. The row name identifies the model added when moving
from \\k-1\\ to \\k\\ candidate models.

## Details

The routine enumerates the available hierarchical log-linear models and
ranks them by BIC. If bootstrap inference is requested, a specified
number of the best-ranked models is retained for bootstrap and jackknife
calculations, after which BCa confidence limits are evaluated for nested
sets of top-BIC models.

This routine implements the bootstrap procedure described by Silverman,
Chan and Vincent (2024). Hierarchical log-linear models are fitted to
the observed data and ordered by increasing BIC.

If `nboot > 0`, multinomial bootstrap samples are generated and the BIC
model-selection procedure is repeated for each bootstrap sample.
Bootstrap model selection is then repeated within nested sets of the
best-ranked models.

With the default settings, all available models are retained for
bootstrap inference for two- and three-list data. For four-list data,
the 20 models with the smallest BIC values are retained, and for five-
and six-list data the 100 models with the smallest BIC values are
retained.

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
<https://doi.org/10.1007/s11222-023-10346-9>.

## Examples

``` r
data(Korea)

# A very small number of bootstrap replications is used here only
# to keep the example quick.
estimate_population_bic(Korea, nboot = 10)
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
#>            0.025      0.1      0.9 0.975
#> [12,23] 134.7273 134.8201 187.6457   201
#> [12,3]  134.7273 134.8201 320.3872   336
#> [13,23] 134.7273 134.8201 320.3872   336
#> [23,1]  134.7273 134.8201 320.3872   336
#> [13,2]  134.7273 134.8201 320.3872   336
#> [1,2,3] 134.7273 134.8201 320.3872   336
#> 
```
