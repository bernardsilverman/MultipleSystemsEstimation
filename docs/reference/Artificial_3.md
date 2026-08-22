# Artificial data set to demonstrate possible instabilities

This is a simple data set based on three lists, which gives examples of
models that fail on one or the other of the criteria tested by
[`check_extended_MLE`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/check_extended_MLE.md).
This is Table 2 in Chan, Silverman and Vincent (2021).

## Usage

``` r
Artificial_3
```

## Format

An object of class `data.frame` with 4 rows and 4 columns.

## Details

If all three two-list effects are included in the fitted model,
[`check_extended_MLE`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/check_extended_MLE.md)
returns status 2: the extended maximum likelihood estimate exists, but
the model matrix is not of full rank, so the parameters are not
identifiable. If the model contains AB either alone or in conjunction
with one of AC and BC, then the extended maximum likelihood estimate
does not exist. If only main effects are considered, or if either or
both of AC and BC, but not AB are included, then the model passes both
tests.

## References

Chan, L., Silverman, B. W., and Vincent, K. (2021). Multiple Systems
Estimation for Sparse Capture Data: Inferential Challenges when there
are Non-Overlapping Lists. *Journal of the American Statistical
Association*, **116(535)**, 1297-1306, Available from
<https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1708748>.
