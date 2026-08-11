# Models BICS, abundance and maxorder.

This routine sorts the models in increasing order according to their
BICs, returns the sorted models with their corresponding BICs and
abundance. The original data as well as the maxorder of the models
considered are returned as well.

## Usage

``` r
assemble_bic(
  xdata,
  maxorder = dim(xdata)[2] - 2,
  checkexist = TRUE,
  removeFRfail = TRUE,
  ...
)
```

## Arguments

- xdata:

  The original data matrix with capture histories and counts.

- maxorder:

  Maximum order of models to be included

- checkexist:

  If TRUE then the Fienberg-Rinaldo condition is checked for each model

- removeFRfail:

  If checkexist is TRUE then models which fail the FR condition are
  removed from the results

- ...:

  Parameters to be fed to `get_hierarchical_models`.

## Value

A list with the following components

- res:

  A matrix with the models considered, their abundance, BIC and their
  order, sorted into increasing order of BIC

- xdata:

  The original data matrix with capture histories and counts.

- maxorder:

  Order parameter that was feed into `get_hierarchical_models`
