# United Kingdom modern-slavery data

Six-list capture-history data from the United Kingdom 2013 strategic
assessment of modern slavery. The lists are local authorities (`LA`),
non-governmental organisations (`NG`), police forces (`PF`), government
organisations (`GO`), the general public (`GP`), and the National Crime
Agency (`NCA`).

## Usage

``` r
UKdat
```

## Format

A data frame with 25 rows and 7 columns. The first six columns are
binary list-membership indicators. Column `count` gives the number of
cases having each observed capture history. Capture histories having
zero count are omitted.

## References

Home Office (2014). Modern Slavery: an application of multiple systems
estimation.
<https://www.gov.uk/government/publications/modern-slavery-an-application-of-multiple-systems-estimation>.
