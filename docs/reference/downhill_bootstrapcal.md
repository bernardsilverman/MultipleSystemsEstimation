# Bootstrap downhill

Construct bootstrap replications and use the downhill fit method to
obtain point estimates of total population sizes from each bootstrap
sample.

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

  original data matrix

- nboot:

  number of bootstrap replicates

- iseed:

  random seed

- checkid:

  If it is TRUE, then `checkident.1` is called and it performs the
  Fienberg-Rinaldo linear program check for the existence of the
  estimates

- verbose:

  If TRUE, return the list of extra output from `downhill_fit`

- maxorder:

  Maximum order of models to be included

## Value

Point estimates of total population sizes from each bootstrap sample.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
<https://doi.org/10.1007/s11222-023-10346-9>.
