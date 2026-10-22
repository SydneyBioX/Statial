# Wrapper to calculate pairwise distance between cell types by image

Calculates the euclidean distance from each cell to the nearest cell of
each type

## Usage

``` r
getDistances(
  cells,
  maxDist = NULL,
  imageID = "imageID",
  spatialCoords = c("x", "y"),
  cellType = "cellType",
  redDimName = "distances",
  distFun = "min",
  nCores = 1
)
```

## Arguments

- cells:

  A dataframe with a cellType column as well as x and y spatial
  coordinates. The dataframe must contain a imageID column and cellID
  (unique cell identifier's) column as well

- maxDist:

  The maximum distance considered.

- imageID:

  The name of the colData column that stores in the image ID.

- spatialCoords:

  The columns that store the spatial coordinates.

- cellType:

  The name of the colData column that stores the cell types.

- redDimName:

  The name of the reduced dimension to store the distances in.

- distFun:

  What distance function to use. Can be min or abundance.

- nCores:

  Number of cores for parallel processing.

## Examples

``` r
data("kerenSCE")

kerenSCE <- getDistances(kerenSCE,
  maxDist = 200
)
```
