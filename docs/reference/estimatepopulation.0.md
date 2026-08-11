# Estimate the total population including the dark figure. If the user wishes to find bootstrap confidence intervals then the routine [`estimatepopulation`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimatepopulation.md) should be used instead.

This routine estimates the total population size, which includes the
dark figure, together with confidence intervals as specified. It also
returns the details of the fitted model. The user can choose whether to
fit main effects only, to fit a particular model containing specified
two-list parameters, or to choose the model using the stepwise approach
described by Chan, Silverman and Vincent (2021).

## Usage

``` r
estimatepopulation.0(
  zdat,
  method = "stepwise",
  quantiles = c(0.025, 0.975),
  mX = NULL,
  pthresh = 0.02
)
```

## Arguments

- zdat:

  Data matrix with \\t+1\\ columns. The first \\t\\ columns, each
  corresponding to a particular list, are 0s and 1s defining the capture
  histories observed. The last column is the count of cases with that
  particular capture history. List names A, B, ... are constructed if
  not supplied. Where a capture history is not explicitly listed, it is
  assumed that it has zero count.

- method:

  If `method = "stepwise"` the stepwise method implemented in
  [`stepwisefit`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/stepwisefit.md)
  is used. If `method = "fixed"` then a specified fixed model is used;
  the model is then given by `mX`. If `method = "main"` then main
  effects only are fitted.

- quantiles:

  Quantiles of interest for confidence intervals.

- mX:

  A \\2 \times k\\ matrix giving the \\k\\ two-list parameters to be
  included in the model if `method = "fixed"`. Each column of `mX`
  contains the numbers of the corresponding pair of lists. If `mX = 0`,
  then all two-list parameters are included. If `mX = NULL`, no two-list
  parameters are included and the main effects model is fitted. If only
  one two-list parameter is to be fitted, it is sufficient to specify it
  as a vector of length 2, e.g `mX=c(1,3)` for the parameter indexed by
  lists 1 and 3. If `method` is equal to `"stepwise"` or `"main"`, then
  `mX` is ignored.

- pthresh:

  Threshold p-value used if `method = "stepwise"`.

## Value

A list with components as below

`estimate` Point estimate and confidence interval estimates
corresponding to specified quantiles.

`MSEfit` The model fitted to the data in the format described in
[`modelfit`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/modelfit.md).

## References

Chan, L., Silverman, B. W., and Vincent, K. (2021). Multiple Systems
Estimation for Sparse Capture Data: Inferential Challenges when there
are Non-Overlapping Lists. *Journal of the American Statistical
Association*, **116(535)**, 1297-1306, Available from
<https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748>.
