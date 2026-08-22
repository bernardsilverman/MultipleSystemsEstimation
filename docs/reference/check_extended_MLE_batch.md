# Carry out the Fienberg-Rinaldo procedure on an array of data vectors and a vector of models

Suppose we have a collection of different data outcomes on the same set
of capture histories and a vector of models. Typically the data outcomes
will be bootstrap replications. This routine finds the unique support
patterns among the data and hence economises the task of finding which
model/data combinations satisfy the Fienberg-Rinaldo condition

## Usage

``` r
check_extended_MLE_batch(x, xcap, zmods)
```

## Arguments

- x:

  a matrix of data observations for a common capture matrix

- xcap:

  the incidence matrix of the capture histories corresponding to the
  rows of x

- zmods:

  a vector of models

## Value

a matrix with rows corresponding to the models and columns to the
columns of x, with elements taking the value TRUE if the FR linear
program for a vector of 0s and 1s with the same zero pattern as the x
data yields a strictly positive value

## References

Silverman, B. W., Chan, L. and Vincent, K., (2024). Bootstrapping
Multiple Systems Estimates to Account for Model Selection *Statistics
and Computing*, **34(44)**, Available from
[\doi{10.1007/s11222-023-10346-9}](NA).
