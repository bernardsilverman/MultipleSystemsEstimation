# Build the model matrix based on particular data, as required to check for identifiability and existence of the maximum likelihood estimate

This routine builds a model matrix as required by the linear program
check
[`check_identifiability`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/check_identifiability.md)
and checks if the matrix is of full rank. In addition, for each
individual list, and for each pair of lists included in the model, it
returns the total count of individuals appearing on the specific list or
lists whether or not in combination with other lists.

## Usage

``` r
buildmodelmatrix(zdat, mX = NULL)
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
  are included. If `mX = NULL`, no such parameters are included and the
  main effects model is fitted.

## Value

A list with components as below

`modmat` The matrix that maps the parameters in the model (excluding any
corresponding to non-overlapping lists) to the log expected value of the
counts of capture histories that do not contain non-overlapping pairs in
the data.

`tvec` A vector indexed by the parameters in the model, excluding those
corresponding to non-overlapping pairs of lists. For each parameter the
vector contains the total count of individuals in all the capture
histories that include both the lists indexed by that parameter.

`rankdef` The column rank deficiency of the matrix `modmat`. If
`rankdef = 0`, the matrix has full column rank.
