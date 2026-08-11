# UK data

Data from the UK 2013 strategic assessment

## Usage

``` r
UKdat
```

## Format

An object of class `data.frame` with 25 rows and 7 columns.

## Details

This is a table of six lists used in the research published by the Home
Office as part of the strategy leading to the Modern Slavery Act 2015.
The data are considered in six lists, labelled as follows: LA–Local
authorities; NG–Non-government organisations; PF–Police forces;
GO–Government organisations; GP–General public; NCA–National Crime
Agency. Each of the first six columns in the data frame corresponds to
one of these lists. Each of the rows of the data frame corresponds to a
possible combination of lists, with value 1 in the relevant column if
the list is in that particular combination. The last column of the data
frame states the number of cases observed in that particular combination
of lists. Combinations of lists for which zero cases are observed are
omitted.

## References

<https://www.gov.uk/government/publications/modern-slavery-an-application-of-multiple-systems-estimation>
