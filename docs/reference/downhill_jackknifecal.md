# Jackknife downhill

Applies the downhill search to the required delete-one data sets and
calculates the BCa acceleration parameter.

## Usage

``` r
downhill_jackknifecal(xdata, checkid = TRUE, maxorder = dim(xdata)[2] - 2)
```

## Arguments

- xdata:

  Capture history data in the standard package format.

- checkid:

  Passed to
  [`downhill_fit()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/downhill_fit.md).

- maxorder:

  Maximum interaction order considered.

## Value

The estimated BCa acceleration parameter.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34**, 44.
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9).
