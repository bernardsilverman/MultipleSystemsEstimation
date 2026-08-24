# Check existence and identifiability of the extended MLE

Checks whether a hierarchical log-linear model satisfies both the
Fienberg–Rinaldo condition for existence of the extended
maximum-likelihood estimate and the rank condition required for
identifiable model parameters.

## Usage

``` r
check_extended_MLE(data, model)
```

## Arguments

- data:

  Capture history data in the standard package format, with one
  indicator column for each list followed by a count column.

- model:

  A character string specifying a hierarchical model, such as
  `"[12,13,23]"`.

## Value

A single integer status code:

- `0`:

  Both conditions are satisfied.

- `1`:

  The Fienberg–Rinaldo condition fails.

- `2`:

  The model parameters are not identifiable.

- `3`:

  Both conditions fail.

## Details

Two distinct conditions are checked.

The Fienberg–Rinaldo condition determines whether the required extended
maximum-likelihood fit exists, accounting for parameters whose
extended-MLE value is minus infinity because the corresponding observed
margins are zero.

The relevant model matrix is also checked for full column rank. Failure
of this condition means that the model parameters are not identifiable.

The two checks are carried out separately, so either one or both may
fail. A model returning a nonzero status should not be fitted using the
package's ordinary maximum-likelihood estimation routines.

## References

Fienberg, S. E. and Rinaldo, A. (2012). Maximum likelihood estimation in
log-linear models. *The Annals of Statistics*, **40**, 996–1023.

Chan, L., Silverman, B. W. and Vincent, K. (2021). Multiple Systems
Estimation for Sparse Capture Data: Inferential Challenges When There
Are Nonoverlapping Lists. *Journal of the American Statistical
Association*, **116**(535), 1297–1306.
[doi:10.1080/01621459.2019.1708748](https://doi.org/10.1080/01621459.2019.1708748).

## Examples

``` r
data(Artificial_3)

check_extended_MLE(
  Artificial_3,
  "[12,13,23]"
)
#> [1] 2

# The result is 2: the existence condition is satisfied,
# but the model parameters are not identifiable.
```
