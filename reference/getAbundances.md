# Wrapper to calculate imhomogenous K function between a cell and surrounding types on each image

Calculate the imhomogenous K function (a measure of cell type abundance)
for each cell to other cell types

## Usage

``` r
getAbundances(
  cells,
  r = 200,
  distFun = "abundance",
  redDimName = "abundances",
  cellType = "cellType",
  imageID = "imageID",
  spatialCoords = c("x", "y"),
  nCores = 1
)
```

## Arguments

- cells:

  A dataframe with a cellType column as well as x and y spatial
  coordinates. The dataframe must contain a imageID column and cellID
  (unique cell identifier's) column as well

- r:

  Radius to include in that calculation of pairwise abundance
  (K-function) between cells (can be a numeric or vector of radii)

- distFun:

  What distance function to use.

- redDimName:

  Name of the reduced dimension to store in sce.

- cellType:

  The name of the column in colData that stores the cell types.

- imageID:

  The name of the column in colData that Stores the image ids.

- spatialCoords:

  The names of the columns in colData that store the spatial
  coordinates.

- nCores:

  Number of cores for parallel processing

## Examples

``` r
library(dplyr)
data("kerenSCE")

singleCellDataCounts <- getAbundances(kerenSCE,
  r = 200,
)
```
