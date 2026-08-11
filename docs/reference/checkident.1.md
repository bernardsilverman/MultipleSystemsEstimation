# The Fienberg-Rinaldo linear program check for the existence of the estimates

This routine performs the Fienberg-Rinaldo linear program check in the
framework of extended maximum likelihood estimates, the parameters
estimates exist if and only if the return value of the check is nonzero

## Usage

``` r
checkident.1(parset, datlist)
```

## Arguments

- parset:

  Either the hierarchical representation of the model, or the vector of
  the corresponding capture histories to the model.

- datlist:

  The output of
  [`ingest_data`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/ingest_data.md)
  on the data set

## Value

The value of the linear program. The parameter estimates within the
extended ML framework exist if and only if this value is nonzero.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2022). Bootstrapping
Multiple Systems Estimates to Account for Model Selection

Fienberg, S. E. and Rinaldo, A. (2012). Maximum likelihood estimation in
log-linear models. *Ann. Statist.* **40**, 996-1023. Supplementary
material: Technical report, Carnegie Mellon University. Available from
<http://www.stat.cmu.edu/~arinaldo/Fienberg_Rinaldo_Supplementary_Material.pdf>.
