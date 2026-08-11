# Comparison of BIC approach and BCa approach

This routine carries out the simulation study as detailed in Section 3.4
of Chan, Silverman and Vincent (2021). If the original data set has low
counts, so that there is a possibility of a simulated data set
containing empty lists, then it may be advisable to use the option
`noninformativelist=TRUE`.

## Usage

``` r
BICandbootstrapsim(
  zdat,
  nsims = 1000,
  nboot = 100,
  pthresh = 0.02,
  iseed = 1234,
  alpha = c(0.025, 0.05, 0.1, 0.16, 0.5, 0.84, 0.9, 0.95, 0.975),
  noninformativelist = FALSE,
  verbose = FALSE,
  ...
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

- nsims:

  Number of simulations to be carried out.

- nboot:

  Number of bootstrap replications for each simulation

- pthresh:

  p-value threshold used in
  [`estimatepopulation.0`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/estimatepopulation.0.md).

- iseed:

  seed for initialization.

- alpha:

  bootstrap quantiles of interests.

- noninformativelist:

  if `noninformativelist=TRUE` then each generated data set in the
  simulation study (including all bootstrap replications) will be
  removed using
  [`tidy_lists`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/tidy_lists.md).

- verbose:

  If `verbose=FALSE`, then the progress of the simulation will not show.
  If `verbose=TRUE`, then the progress of the simulation will be shown.

- ...:

  other arguments.

## Value

A list with components as below

`popest` Total population point estimate from the original data using
`estimatepopulation.0` with default threshold.

`BICmodels` The best model chosen by the BIC at each simulation.

`BICvals` Point estimates of the total population and standard error of
the best model chosen by the BIC at each simulation.

`simreps` Counts associated to each capture history at each simulation.

`modelmat` A full capture history matrix excluding the row corresponding
to the dark figure.

`popestsim` Total population estimate given by the BCa method in each
simulation.

`BCaquantiles` bootstrap confidence intervals given by the BCa method.

`BICconf` confidence interval given by the BIC method.

## Details

This function requires the optional package Rcapture. If it is installed
but not attached, the user will be prompted to load it when the function
is called. If it is not installed, the function will display
installation instructions and return without running the simulation.

## References

Chan, L., Silverman, B. W., and Vincent, K. (2021). Multiple Systems
Estimation for Sparse Capture Data: Inferential Challenges when there
are Non-Overlapping Lists. *Journal of the American Statistical
Association*, **116(535)**, 1297-1306, Available from
<https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748>.

DiCiccio, T. J. and Efron, B. (1996). Bootstrap Confidence Intervals.
*Statistical Science*, **40(3)**, 189-228.

Rivest, L-P. and Baillargeon, S. (2014) Rcapture. CRAN package.
Available from Available from
<https://CRAN.R-project.org/package=Rcapture>.
