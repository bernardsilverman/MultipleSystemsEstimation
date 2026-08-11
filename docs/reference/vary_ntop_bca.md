# BCa inference for varying numbers of top BIC models

Investigates the effect on bootstrap inference of restricting model
selection in the bootstrap and jackknife calculations to a subset of the
candidate models.

The population-size point estimate is obtained by minimizing BIC over
the full set of candidate models. For the bootstrap and jackknife
calculations, however, model selection may be restricted to the `ntop`
models having the smallest BIC values for the original data. This can
greatly reduce the computational burden.

This function calculates BCa inference for all values of `ntop` up to a
specified maximum, allowing the effect of this restriction on bootstrap
inference to be investigated. The maximum may be set to `Inf`, in which
case all candidate models are considered.

The models are ranked using the original-data BIC values. For each
bootstrap and jackknife replication, and for each value of `ntop` up to
the specified maximum, the best-fitting model among the first `ntop`
models in this ranking is selected. The calculations for the different
values of `ntop` are obtained from the same set of model fits.

With the default `degree = 1`, this function implements the algorithm
described in Section 2.5 of Silverman, Chan and Vincent (2024), and can
be used to reproduce the calculations underlying Figures 1 and 2 of that
paper. An optional neighbourhood-based ranking, corresponding to
`degree = 2` and considered in Section 5.1, is also available; see
Details.

## Usage

``` r
vary_ntop_bca(
  zdat,
  maxorder,
  ntopmax = 50,
  degree = 1,
  nboot = 1000,
  iseed = 1234,
  alpha = c(0.025, 0.1, 0.9, 0.975)
)
```

## Arguments

- zdat:

  The capture data, in the standard format used by
  MultipleSystemsEstimation.

- maxorder:

  The maximum order of interaction allowed in the hierarchical loglinear
  models considered. Must be at least 2.

- ntopmax:

  The largest value of `ntop` to be considered. Inference is calculated
  for every integer value of `ntop` from 1 to `ntopmax`. If
  `ntopmax = Inf`, all candidate models are considered and inference is
  calculated for every value of `ntop` up to the total number of
  candidate models.

- degree:

  The degree of the model ranking. Must be either 1 or 2. The default,
  `degree = 1`, orders models by their BIC values for the original data.
  The option `degree = 2` uses the alternative neighbourhood-based
  ranking investigated in Section 5.1 of Silverman, Chan and Vincent
  (2024).

- nboot:

  The number of bootstrap replications.

- iseed:

  The random-number seed used for the bootstrap.

- alpha:

  The probabilities at which BCa confidence limits are required.

## Value

A list with two components:

- `estimate`:

  The population-size point estimate obtained by minimizing BIC over the
  full set of candidate models.

- `inference`:

  A data frame with one row for each value of `ntop` considered. The
  first column gives `ntop`; the remaining columns give the requested
  BCa confidence limits. If `ntopmax = Inf`, rows are returned for every
  value of `ntop` up to the total number of candidate models.

## Details

For each bootstrap and jackknife replication, the population-size
estimate and BIC are calculated for every model up to `ntopmax` in the
chosen ranking. Results for smaller values of `ntop` are then obtained
from these fitted models without repeating the model fits.

If `ntopmax = Inf`, all candidate models are fitted for each bootstrap
and jackknife replication. In this case the ordering has no effect on
the resulting inference, because all candidate models are considered.
The exhaustive calculation may be computationally expensive when the
number of candidate models is large.

With `degree = 2`, the degree-2 BIC rank of a model is the smallest
ordinary BIC rank among its 1-neighbours, including the model itself,
where two models are 1-neighbours if the symmetric difference between
their sets of parameters has size at most 1. Ties in degree-2 rank are
broken by ordinary BIC rank.

The `degree = 2` option is included mainly to reproduce and investigate
the alternative ordering considered in Section 5.1 of Silverman, Chan
and Vincent (2024). In that paper the neighbourhood-based ordering did
not improve on the ordinary BIC ordering in a useful way.

The results are intended to inform the choice of `ntop` used in the
user-facing BIC bootstrap routines, in particular the default and
guidance for
[`estimate_population_bic`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimate_population_bic.md).
By comparing inference obtained with restricted values of `ntop` against
inference obtained using larger values, the calculation provides
empirical evidence about how much computational economy can be gained
without materially changing the resulting confidence intervals.

## References

- Silverman, B. W., Chan, L. and Vincent, K. (2024). Bootstrapping
  multiple systems estimates to account for model selection. *Statistics
  and Computing*, **34**, 44.
  [doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9)
  .

- Efron, B. and Tibshirani, R. (1986). Bootstrap methods for standard
  errors, confidence intervals, and other measures of statistical
  accuracy. *Statistical Science*, **1**, 54–75.
