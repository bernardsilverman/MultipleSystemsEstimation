# Subset matrix

The idea of this routine is to reduce either or both of `ntopmodels` and
`maxorder` without the need to recalculate any actual model fits.

## Usage

``` r
subsetmat(z, ntopmodels = Inf, maxorder = Inf)
```

## Arguments

- z:

  output from `assemble_bic` or a list output from applying
  `assemble_bic`, `jackknifecal` and `bootstrapcal`.

- ntopmodels:

  number of top models. If (taking into account any change in the
  maximum order of models) there are fewer than `ntopmodels` in the data
  supplied, then it will be reduced to that value. If it is not
  specified then there will be no reduction in the number.

- maxorder:

  the maximum order of the models to be considered. If not specified, it
  will be set to the corresponding value in the input data.

## Value

A list with the following components

- res:

  a matrix containing models being considered, abundance, BIC and their
  ordered after being subsetted by maxorder and ntopmodels

- xdata:

  Original data matrix with counts and capture histories

- maxorder:

  The maximum order of models considered after subsetting

- jackabund:

  Jackknife abundance matrix, subsetted by maxorder and ntopmodels

- jackbic:

  Jackknife BIC matrix, subsetted by maxorder and ntopmodels

- countsobserved:

  Capture counts in the same order as the columns of `jackabund` and
  `jackbic`

- bootabund:

  Bootstrap abundance matrix, subsetted by maxorder and ntopmodels

- bootbic:

  Bootstrap BIC matrix, subsetted by maxorder and ntopmodels

If the input only has the output from `assemble_bic`, the last five
items of the list do not appear.

## Details

The routine subsets the results matrix as part of the output given by
`assemble_bic` based on specified parameters `ntopmodels` and
`maxorder`. It returns the subsetted matrix, original data matrix with
capture histories and counts and the new actual value of `maxorder`
(reducing it from the input value if necessary or if the default input
value of \\\infty\\ is used).
