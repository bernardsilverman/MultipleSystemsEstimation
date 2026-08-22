# Bootstrap abundance and bic

This routine takes the output from `assemble_bic` or `subsetmat` and
returns bootstrap abundance matrix and BIC matrix. This version makes
use of `check_extended_MLE_batch`.

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

  Results from `assemble_bic` or `subsetmat`

- nboot:

  The number of bootstrap replications.

- iseed:

  Integer seed to allow for replicability.

- checkexist:

  If `checkexist=TRUE`, check for existence, else it does not check for
  existence.

- saveinterval:

  If this is set to a finite value, the output list `z` will be saved
  every time the number of replications is a multiple of it. A message
  will be printed every time it is saved.

- savefile:

  The file to which the output will be saved if `saveinterval` is set to
  a finite value.

## Value

The original input list `z` with the additional components

- bootabund:

  Bootstrap abundance matrix

- bootbic:

  Bootstrap BIC matrix

If there are already components with these names they will be
overwritten.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
[\doi{10.1007/s11222-023-10346-9}](NA).
