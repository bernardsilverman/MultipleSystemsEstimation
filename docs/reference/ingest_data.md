# Preliminary processing of a data matrix

Perform various preprocessing tasks on the data

## Usage

``` r
ingest_data(xdat)
```

## Arguments

- xdat:

  Data matrix of the usual kind

## Value

A list with the following elements

- nobs:

  Numbers of observations indexed by encoded histories

- nstar:

  For each capture history, total number of observations for that
  capture history and all its descendants

- nlists:

  Total number of lists

- listnames:

  Names of the lists, constructed to be A, B, ... if necessary

- data:

  The input data matrix

- notestimable:

  A vector indicating which parameters are not estimable, because they
  are strict descendants of parameters which would be already estimated
  to be \\-\infty\\ if they are included in the model

- masterdesign:

  The inclusion matrix as constructed by
  [`make_master_design`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/make_master_design.md)

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
[\doi{10.1007/s11222-023-10346-9}](NA).
