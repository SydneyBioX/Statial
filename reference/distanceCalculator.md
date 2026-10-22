# Calculate pairwise distance between cell types

Calculates the euclidean distance from each cell to the nearest cell of
each type for a single image

## Usage

``` r
distanceCalculator(data, maxDist = 200, distFun = "min")
```

## Arguments

- data:

  the single cell data of interest

- maxDist:

  Maximum distance between pairs of points to be counted as close pairs.

- distFun:

  How to merge duplicate entries.
