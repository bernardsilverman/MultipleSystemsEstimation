# Jackknife downhill

It uses the downhill approach to calculate the jackknife abundance and
returns the estimated acceleration factor

## Usage

``` r
downhill_jackknifecal(xdata, checkid = TRUE, maxorder = dim(xdata)[2] - 2)
```

## Arguments

- xdata:

  original data matrix

- checkid:

  If it is TRUE, then `checkident.1` is called and it performs the
  Fienberg-Rinaldo linear program check for the existence of the
  estimates

- maxorder:

  Maximum order of models to be included

## Value

the estimated acceleration factor

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
<https://doi.org/10.1007/s11222-023-10346-9>.
