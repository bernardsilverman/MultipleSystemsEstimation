# Bootstrap downhill

Generates multinomial bootstrap samples and applies
[`downhill_fit()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/downhill_fit.md)
to each one.

## Usage

``` r
downhill_bootstrapcal(
  xdata,
  nboot = 1000,
  iseed = 1234,
  checkid = TRUE,
  verbose = FALSE,
  maxorder = dim(xdata)[2] - 2
)
```

## Arguments

- xdata:

  Capture history data in the standard package format.

- nboot:

  Number of bootstrap replications.

- iseed:

  Integer random-number seed.

- checkid:

  Passed to
  [`downhill_fit()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/downhill_fit.md).

- verbose:

  Passed to
  [`downhill_fit()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/downhill_fit.md).

- maxorder:

  Maximum interaction order considered.

## Value

If `verbose = FALSE`, a numeric vector of bootstrap population
estimates. If `verbose = TRUE`, the detailed results returned by
[`downhill_fit()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/downhill_fit.md)
for each replication.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34**, 44.
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9).
