# Calculate the level of marker contamination of each cell

Calculates contamination scores using a random forest classification

## Usage

``` r
calcContamination(
  cells,
  markers = NULL,
  num.trees = 100,
  verbose = FALSE,
  missingReplacement = 0,
  assay = "intensities",
  cellType = "cellType",
  redDimName = "contaminations"
)
```

## Arguments

- cells:

  A SingleCellExperiment or SpatialExperiment with a cellType column as
  well as marker intensity information corresponding to each cell.

- markers:

  A vector of markers that proxy a cell's state. If NULL, all markers
  will be used.

- num.trees:

  Number of trees to be used in the random forest classifier

- verbose:

  A logical indicating whether information about the final random forest
  model should be outputted.

- missingReplacement:

  A default value to replace missing marker intensities for
  classification.

- assay:

  The assay in the SingleCellExperiment object that contains the desired
  marker expressions.

- cellType:

  The name of the column in colData that stores the cell types.

- redDimName:

  The redDimName to store the output in the sce.

## Examples

``` r
data("kerenSCE")

singleCellDataDistancesContam <- calcContamination(
  kerenSCE
)
```
