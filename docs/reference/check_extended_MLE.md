# Check existence and identifiability of the extended MLE for a given model and dataset

Checks whether a specified hierarchical log-linear model satisfies both
the Fienberg–Rinaldo condition for maximum-likelihood estimation and the
rank condition required for identifiable model parameters.

## Usage

``` r
check_extended_MLE(data, model)
```

## Arguments

- data:

  A data frame or matrix containing the capture histories and their
  observed counts. The first columns are binary indicators for the
  capture lists and the final column contains the count for each capture
  history. Missing capture histories are treated as having count zero.

- model:

  The model to be checked. This may be either a character string giving
  the model in hierarchical notation, such as `"[12,13,23]"`, or a
  numeric vector of encoded model parameters, such as that returned by
  [`convert_from_hierarchy`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/convert_from_hierarchy.md).

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

A model returning a nonzero status should not be fitted using the
package's ordinary maximum-likelihood estimation routines.

## Details

Two distinct conditions are checked.

First, the Fienberg–Rinaldo condition determines whether the required
maximum-likelihood fit exists, after accounting for parameters forced to
the boundary by zero observed margins.

Second, the relevant model matrix is checked for full column rank.
Failure of this condition means that the model parameters are not
identifiable.

The model may be supplied in hierarchical notation or as its
corresponding vector of encoded parameters. Hierarchical notation
specifies the generators of the model; all lower-order terms required by
hierarchy are included automatically. See
[`convert_from_hierarchy`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/convert_from_hierarchy.md).

The two checks are carried out independently, so it is possible for
either one or both to fail.

## References

Fienberg, S. E. and Rinaldo, A. (2012). Maximum likelihood estimation in
log-linear models. *The Annals of Statistics*, **40**, 996–1023.

Chan, L., Silverman, B. W. and Vincent, K. (2021). Multiple systems
estimation for sparse capture data: Inferential challenges when there
are nonoverlapping lists. *Journal of the American Statistical
Association*, **116**, 1297–1310.

## Examples

``` r
data(Artificial_3)

# Specify the model in hierarchical notation.
check_extended_MLE(
  Artificial_3,
  "[12,13,23]"
)
#> [1] 2

# Equivalently, use the encoded parameter representation.
encoded_model <- convert_from_hierarchy("[12,13,23]")
#> Error in convert_from_hierarchy("[12,13,23]"): could not find function "convert_from_hierarchy"

check_extended_MLE(
  Artificial_3,
  encoded_model
)
#> Error: object 'encoded_model' not found

# Both calls return 2: the existence condition is satisfied,
# but the model parameters are not identifiable.
```
