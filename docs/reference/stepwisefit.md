# Stepwise fit using Poisson p-values.

Starting with a model with main effects only, two-list parameters are
added one by one. At each stage the parameter with the lowest p-value is
added, provided that p-value is lower than `pthresh`, and provided that
the resulting model does not fail either of the tests in
[`check_identifiability`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/check_identifiability.md).

## Usage

``` r
stepwisefit(zdat, pthresh = 0.02)
```

## Arguments

- zdat:

  Data matrix with \\t+1\\ columns. The first \\t\\ columns, each
  corresponding to a particular list, are 0s and 1s defining the capture
  histories observed. The last column is the count of cases with that
  particular capture history. List names A, B, ... are constructed if
  not supplied. Where a capture history is not explicitly listed, it is
  assumed that it has zero count.

- pthresh:

  this is the threshold below which the p-value of the newly added
  parameter needs to be in order to be included in the model. If
  `pthresh = 0` then the model with main effects only is returned.

## Value

A list with components as below

`fit` Details of the fit of the specified model as output by `glm`. The
Akaike information criterion is adjusted to take account of the number
of parameters corresponding to non-overlapping pairs.

`emptyoverlaps` Matrix with two rows, giving the list pairs within the
model for which no cases are observed in common. Each column gives the
indices of a pair of lists, with the names of the lists in the column
name.

`poisspempty` the Poisson p-values of the parameters corresponding to
non-overlapping pairs.

## Details

For each candidate two-list parameter for possible addition to the
model, the p-value is calculated as follows. The total of cases
occurring on both lists indexed by the parameter (regardless of whether
or not they are on any other lists) is calculated. On the null
hypothesis that the effect is not included in the model, this statistic
has a Poisson distribution whose mean depends on the parameters within
the model. The one-sided Poisson p-value of the observed statistic is
calculated.
