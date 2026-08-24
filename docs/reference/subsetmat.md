# Restrict a set of fitted hierarchical models

Restricts previously calculated model fits by maximum interaction order
and original-data BIC rank, without repeating any fits. Bootstrap and
jackknife matrices already present in the input are subsetted in
parallel.

## Usage

``` r
subsetmat(z, ntopmodels = Inf, maxorder = Inf)
```

## Arguments

- z:

  A result from
  [`assemble_bic()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/assemble_bic.md),
  optionally augmented by
  [`bootstrapcal()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/bootstrapcal.md)
  and
  [`jackknifecal()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/jackknifecal.md).

- ntopmodels:

  Maximum number of models to retain after applying the
  interaction-order restriction. The default `Inf` retains all available
  models.

- maxorder:

  Maximum interaction order to retain. The default `Inf` imposes no
  additional restriction.

## Value

The input list `z`, subsetted to the retained models. Its possible
components are:

- `res`:

  The original-data model results.

- `xdata`:

  The original capture history data.

- `maxorder`:

  The largest retained interaction order.

- `jackabund`, `jackbic`:

  Jackknife population estimates and BIC values, if present in the
  input.

- `countsobserved`:

  Capture history counts, if present in the input.

- `bootabund`, `bootbic`:

  Bootstrap population estimates and BIC values, if present in the
  input.
