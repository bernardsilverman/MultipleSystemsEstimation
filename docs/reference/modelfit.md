# Fit a specified model to multiple systems estimation data

This routine fits a specified model to multiple systems estimation data,
taking account of the possibility of empty overlaps between pairs of
observed lists.

## Usage

``` r
modelfit(zdat, mX = NULL, check = TRUE)
```

## Arguments

- zdat:

  Data matrix with \\t+1\\ columns. The first \\t\\ columns, each
  corresponding to a particular list, are 0s and 1s defining the capture
  histories observed. The last column is the count of cases with that
  particular capture history. List names A, B, ... are constructed if
  not supplied. Where a capture history is not explicitly listed, it is
  assumed that it has zero count.

- mX:

  A \\2 \times k\\ matrix giving the \\k\\ two-list parameters to be
  included in the model. Each column of `mX` contains the numbers of the
  corresponding pair of lists. If `mX = 0`, then all two-list parameters
  are included. If `mX = NULL`, no two-list parameters are included and
  the main effects model is fitted. If only one two-list parameter is to
  be fitted, it may be specified as a vector of length 2, e.g
  `mX=c(1,3)` for the parameter corresponding to lists 1 and 3.

- check:

  If `check = TRUE` check first of all if the maximum likelihood
  estimate exists and is identifiable, using the routine
  [`check_identifiability`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/check_identifiability.md).
  If either condition fails, print an appropriate error message and
  return the error code.

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
