# Artificial three-list data

An artificial three-list data set used to demonstrate failures of the
conditions tested by
[`check_extended_MLE`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/check_extended_MLE.md).
The data appear in Table 2 of Chan, Silverman and Vincent (2021).

## Usage

``` r
Artificial_3
```

## Format

A data frame with 4 rows and 4 columns. Columns `A`, `B` and `C` are
binary list-membership indicators. Column `n` gives the number of cases
having each observed capture history.

## Details

If all three two-list interactions are included,
[`check_extended_MLE()`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/check_extended_MLE.md)
returns status 2: the extended MLE exists, but the model matrix is not
of full rank and the parameters are not identifiable. If the model
contains AB, either alone or together with AC or BC, the extended MLE
does not exist. The main-effects model and models containing either or
both of AC and BC, but not AB, pass both checks.

## References

Chan, L., Silverman, B. W. and Vincent, K. (2021). Multiple Systems
Estimation for Sparse Capture Data: Inferential Challenges When There
Are Nonoverlapping Lists. *Journal of the American Statistical
Association*, **116**(535), 1297–1306.
[doi:10.1080/01621459.2019.1708748](https://doi.org/10.1080/01621459.2019.1708748).
