# Build model for multiple systems estimation

For multiple systems estimation model corresponding to a specified set
of two-list effects, set up the GLM model formula and data matrix.

## Usage

``` r
buildmodel(zdat, mX)
```

## Arguments

- zdat:

  Data matrix with \\t+1\\ columns. The first \\t\\ columns, each
  corresponding to a particular list, are 0s and 1s defining the capture
  histories observed. The last column is the count of cases with that
  particular capture history. List names A, B, ... are constructed if
  not supplied. Where a capture history is not explicitly listed, it is
  assumed that it has zero counts.

- mX:

  A \\2 \times k\\ matrix giving the \\k\\ two-list effects to be
  included in the model. Each column of `mX` contains the numbers of the
  corresponding pair of lists. If `mX = 0`, then all two-list effects
  are included. If `mX = NULL`, no such effects are included and the
  main effects model is fitted.

## Value

A list with components as below.

`datamatrix` A matrix with all possible capture histories, other than
those equal to or containing non-overlapping pairs indexed by parameters
that are within the model specified by `mX`. A non-overlapping pair is a
pair of lists \\(i,j)\\ such that no case is observed in both lists,
regardless of whether it is present on any other lists. If \\(i,j)\\ is
within the model specified by `mX`, all capture histories containing
both \\i\\ and \\j\\ are then excluded.

`modelform` The model formula suitable to be called by the Generalized
Linear Model function `glm`. Model terms corresponding to
non-overlapping pairs are not included, because they are handled by
removing appropriate rows from the data matrix supplied to glm. The list
of non-overlapping pairs are provided in `emptyoverlaps`. See Chan,
Silverman and Vincent (2021) for details.

`emptyoverlaps` A matrix with two rows, whose columns give the indices
of non-overlapping pairs of lists where the parameter indexed by the
pair is within the specified model. The column names give the names of
the lists corresponding to each pair.

## References

Chan, L., Silverman, B. W., and Vincent, K. (2021). Multiple Systems
Estimation for Sparse Capture Data: Inferential Challenges when there
are Non-Overlapping Lists. *Journal of the American Statistical
Association*, **116(535)**, 1297-1306, Available from
<https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748>.
