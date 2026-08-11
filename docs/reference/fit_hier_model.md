# Fit a hierarchical model taking account of possible sparsity

Fit a hierarchical model taking account of possible sparsity

## Usage

``` r
fit_hier_model(xdatin, hiermod, bicRcap = TRUE, checkid = FALSE)
```

## Arguments

- xdatin:

  data obtained using `ingest_data`

- hiermod:

  hierarchical model to fit

- bicRcap:

  if TRUE then use the Rcapture convention that the BIC sample size is
  the number of cases observed. Otherwise use the number of cells in the
  Poisson log linear model.

- checkid:

  if TRUE then `checkident.1` is called inside the routine

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
<https://doi.org/10.1007/s11222-023-10346-9>.
