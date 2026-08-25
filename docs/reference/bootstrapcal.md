# Bootstrap abundance and BIC values

Generates multinomial bootstrap samples and fits every retained
hierarchical model to each sample. When requested, model-data
combinations are screened using
[`check_extended_MLE_batch()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/check_extended_MLE_batch.md)
before fitting.

## Usage

``` r
bootstrapcal(
  z,
  nboot = 1000,
  iseed = 1234,
  checkexist = TRUE,
  saveinterval = Inf,
  savefile = "bootout.Rdata"
)
```

## Arguments

- z:

  A result from
  [`assemble_bic()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/assemble_bic.md)
  or
  [`subsetmat()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/subsetmat.md).

- nboot:

  Number of bootstrap replications.

- iseed:

  Integer random-number seed.

- checkexist:

  If `TRUE`, check identifiability and existence of the extended MLE for
  each model and support pattern before fitting.

- saveinterval:

  If finite, save the accumulating result whenever the replication
  number is a multiple of this value.

- savefile:

  File used when `saveinterval` is finite.

## Value

The input list `z`, with the following components added or replaced:

- `countsobserved`:

  Counts for the complete set of observable capture histories.

- `bootreplications`:

  A matrix whose columns contain the bootstrap counts.

- `bootabund`:

  A matrix of population estimates, with models in rows and bootstrap
  replications in columns.

- `bootbic`:

  A corresponding matrix of BIC values.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection. *Statistics
and Computing*, **34**, 44.
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9)
.
