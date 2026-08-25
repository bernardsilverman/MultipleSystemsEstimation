# Check extended-MLE conditions for multiple data vectors and models

Suppose we have a vector of models and a collection of different data
outcomes on the same set of capture histories. Typically, these will be
bootstrap replications. This routine carries out the extended-MLE checks
for every combination of data outcome and model. It economises the task
of determining which model/data combinations satisfy both the
identifiability and Fienberg–Rinaldo conditions by first finding the
unique support patterns among the data outcomes.

## Usage

``` r
check_extended_MLE_batch(x, xcap, zmods)
```

## Arguments

- x:

  Numeric matrix whose columns contain count vectors for a common set of
  capture histories.

- xcap:

  Binary matrix defining the capture histories corresponding to the rows
  of `x`.

- zmods:

  Character vector of hierarchical-model strings.

## Value

A logical matrix with models in rows and data vectors in columns. An
element is `TRUE` when both extended-MLE conditions are satisfied for
that model and support pattern.

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34**, 44,
[doi:10.1007/s11222-023-10346-9](https://doi.org/10.1007/s11222-023-10346-9)
.
