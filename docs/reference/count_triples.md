# Count number of triples of overlapping lists

The routine counts the number of subsets of size three of lists such
that every pair of lists in the triple overlaps. If the number is zero,
then the model with all two-list effects is unidentifiable.

## Usage

``` r
count_triples(zdat)
```

## Arguments

- zdat:

  Data matrix with \\t+1\\ columns. The first \\t\\ columns, each
  corresponding to a particular list, are 0s and 1s defining the capture
  histories observed. The last column is the count of cases with that
  particular capture history. List names A, B, ... are constructed if
  not supplied. Where a capture history is not explicitly listed, it is
  assumed that it has zero count.

## Value

a count of subsets of size three of lists such that every pair of lists
in the triple overlaps.
