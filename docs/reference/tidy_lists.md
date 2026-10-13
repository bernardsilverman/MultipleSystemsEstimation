# Consolidate capture history data

Combines duplicate capture histories by summing their counts. It can
also add all observable histories having zero count and remove
non-informative capture lists.

## Usage

``` r
tidy_lists(zdat, includezerocounts = FALSE, remove_noninformative = FALSE)
```

## Arguments

- zdat:

  A matrix or data frame with \\t+1\\ columns. The first \\t\\ columns
  are binary list-membership indicators and the final column contains
  the capture history counts. Histories not included explicitly are
  assumed to have zero count.

- includezerocounts:

  If `TRUE`, include every observable capture history, including those
  with zero count. The all-zero history representing the unobserved
  population is never included. If `FALSE`, return only positive-count
  histories.

- remove_noninformative:

  If `TRUE`, remove lists containing all or none of the observed cases
  and remove duplicate list columns, retaining the first copy.

## Value

A data frame with one row for each retained capture history and one
column for each retained list, followed by the count column.

## Examples

``` r
# Two rows have the same capture history: A only.
zdat <- data.frame(
  A = c(1, 1, 0, 1),
  B = c(0, 0, 1, 1),
  C = c(0, 0, 0, 1),
  count = c(3, 5, 4, 2)
)
zdat
#>   A B C count
#> 1 1 0 0     3
#> 2 1 0 0     5
#> 3 0 1 0     4
#> 4 1 1 1     2

# Consolidate rows with the same capture history by summing their counts.
# The two A-only rows become one row with count 3 + 5 = 8.
tidy_lists(zdat)
#>   A B C count
#> 1 1 0 0     8
#> 2 0 1 0     4
#> 3 1 1 1     2

# Include all seven observable histories, filling missing counts with zero.
# The unobserved history (0, 0, 0) is excluded.
tidy_lists(zdat, includezerocounts = TRUE)
#>   A B C count
#> 1 1 0 0     8
#> 2 0 1 0     4
#> 3 0 0 1     0
#> 4 1 1 0     0
#> 5 1 0 1     0
#> 6 0 1 1     0
#> 7 1 1 1     2

# Add a duplicate list, an empty list and a list containing every case.
# Keep the count column last.
redundant <- data.frame(
  zdat[c("A", "B", "C")],
  A_copy = zdat$A,
  empty = 0,
  everyone = 1,
  count = zdat$count
)
redundant
#>   A B C A_copy empty everyone count
#> 1 1 0 0      1     0        1     3
#> 2 1 0 0      1     0        1     5
#> 3 0 1 0      0     0        1     4
#> 4 1 1 1      1     0        1     2

# Remove the three redundant lists and consolidate repeated histories.
tidy_lists(redundant, remove_noninformative = TRUE)
#>   A B C count
#> 1 1 0 0     8
#> 2 0 1 0     4
#> 3 1 1 1     2
```
