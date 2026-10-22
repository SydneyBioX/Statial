# Evaluation of pairwise cell relationships, conditional on a 3rd population.

Kontextual identifies the relationship between two cell types which are
conditional on the spatial behaviour of a 3rd cell population, for a
particular radius (r).

## Usage

``` r
Kontextual(
  cells,
  r,
  parentDf = NULL,
  from = NULL,
  to = NULL,
  parent = NULL,
  image = NULL,
  inhom = FALSE,
  edgeCorrect = TRUE,
  window = "convex",
  window.length = NA,
  includeOriginal = TRUE,
  spatialCoords = c("x", "y"),
  cellType = "cellType",
  imageID = "imageID",
  cores = 1
)
```

## Arguments

- cells:

  A SingleCellExperiment, SpatialExperiment or a list of data.frames
  containing columns specifying the imageID, cellType, and x and y
  spatial coordinates.

- r:

  Radii to evaluated pairwise relationships between from and to cells.

- parentDf:

  A data frame from
  [`parentCombinations`](https://sydneybiox.github.io/Statial/reference/parentCombinations.md)

- from:

  The first cell type to be evaluated in the pairwise relationship.

- to:

  The second cell type to be evaluated in the pairwise relationship.

- parent:

  The parent population of the from cell type (must include from cell
  type).

- image:

  A vector of images to subset the results to. If NULL we default to all
  images.

- inhom:

  A logical value indicating whether to account for inhomogeneity.

- edgeCorrect:

  A logical value indicating whether to perform edge correction.

- window:

  Type of window for data, either \`square\`, \`convex\` or \`concave\`,
  passed into
  [`makeWindow`](https://sydneybiox.github.io/Statial/reference/makeWindow.md)

- window.length:

  A tuning parameter for controlling the level of concavity when
  estimating concave windows. Passed into
  [`makeWindow`](https://sydneybiox.github.io/Statial/reference/makeWindow.md)

- includeOriginal:

  A logical value to return the original L function values along with
  the kontextual values.

- spatialCoords:

  The columns which contain the x and y spatial coordinates.

- cellType:

  The column which contains the cell types.

- imageID:

  The column which contains image identifiers.

- cores:

  Number of cores for parallel processing.

## Value

A kontextualResult object

## Examples

``` r
# Load data
data("kerenSCE")


CD4_Kontextual <- Kontextual(
  cells = kerenSCE,
  r = 50,
  from = "Macrophages",
  to = "Keratin_Tumour",
  parent = c("Macrophages", "CD4_Cell"),
  image = "6"
)


head(CD4_Kontextual)
#>   imageID                        test  original kontextual  r inhomL
#> 1       6 Macrophages__Keratin_Tumour -25.81863  -33.00758 50  FALSE
```
