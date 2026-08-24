# Five-list version of the United Kingdom data

A five-list version of
[`UKdat`](https://bernardsilverman.github.io/MultipleSystemsEstimation/reference/UKdat.md),
constructed by combining the police-force and National Crime Agency
lists into a single list, `PFNCA`.

## Usage

``` r
UKdat_5
```

## Format

A data frame with 18 rows and 6 columns. Columns `LA`, `NG`, `PFNCA`,
`GO` and `GP` are binary list-membership indicators. Column `count`
gives the number of cases having each observed capture history.
