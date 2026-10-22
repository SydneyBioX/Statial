# First layer wrapper function to build linear models measuring state changes

Builds linear models measuring marker based state changes in a cell type
based of the proximity or abundance of another cell type. The function
provides the option to build robust and mixed linear model variants

## Usage

``` r
calcStateChanges(
  cells,
  marker = NULL,
  from = NULL,
  to = NULL,
  image = NULL,
  type = "distances",
  assay = 1,
  cellType = "cellType",
  imageID = "imageID",
  contamination = NULL,
  test = "g",
  minCells = 20,
  verbose = FALSE,
  timeout = 10,
  nCores = 1
)
```

## Arguments

- cells:

  A dataframe with a imageID, cellType, and marker intensity column
  along with covariates (e.g. distance or abundance of the nearest cell
  type) to model cell state changes

- marker:

  A vector of markers that proxy a cell's state. If NULL, all markers
  will be used.

- from:

  A vector of cell types to use as the primary cells. If NULL, all cell
  types will be used.

- to:

  A vector of cell types to use as the interacting cells. If NULL, all
  cell types will be used.

- image:

  A vector of images to filter to. If null all images will be used.

- type:

  What type of state change. This value should be in reduced dimensions.

- assay:

  The assay in the SingleCellExperiment object that contains the marker
  expressions.

- cellType:

  The column in colData that stores the cell types.

- imageID:

  The column in colData that stores the image ids.

- contamination:

  If TRUE, use the contamination scores that have previously been
  calculate. Otherwise a name of which reduced dimension contains the
  scores.

- test:

  The type of test to perform. By default this will assume the data is
  Gaussian. A value of "nb" will use a negative binomial to model the
  expression.

- minCells:

  The minimum number of cells required to fit a model.

- verbose:

  A logical indicating if messages should be printed

- timeout:

  A maximum time allowed to build each model. Setting this may be
  important when building rlm mixed linear models

- nCores:

  Number of cores for parallel processing

## Examples

``` r
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
data("kerenSCE")

kerenSCE <- kerenSCE[, kerenSCE$imageID %in% c(5, 6)]

kerenSCE <- getDistances(kerenSCE,
  maxDist = 200,
)

imageModels <- calcStateChanges(
  cells = kerenSCE,
  from = "Macrophages",
  to = "Tumour"
)
```
