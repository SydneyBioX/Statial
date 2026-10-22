# Evaluation of Kontextual over a range of radii.

This function obtains \`Kontextual\` values over a range of radii,
standard deviations for each value can be obtained using permutation for
significance testing. To obtain estimates for standard deviations
specify \`se = TRUE\`.

## Usage

``` r
kontextCurve(
  cells,
  from,
  to,
  parent,
  image = NULL,
  rs = seq(10, 100, 10),
  inhom = FALSE,
  edge = TRUE,
  se = FALSE,
  nSim = 20,
  cores = 1,
  imageID = "imageID",
  cellType = "cellType",
  ...
)
```

## Arguments

- cells:

  A single image from a SingleCellExperiment object

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

- rs:

  A vector of radii to evaluate kontextual over.

- inhom:

  A logical value indicating whether to perform an inhomogeneous L
  function.

- edge:

  A logical value indicating whether to perform edge correction.

- se:

  A logical value to indicate if the standard deviation of kontextual
  should be calculated to construct error bars.

- nSim:

  Number of randomisations to perform using
  [`relabelKontextual`](https://sydneybiox.github.io/Statial/reference/relabelKontextual.md),
  which will be used to calculated the SE.

- cores:

  Number of cores for parallel processing.

- imageID:

  The column in colData that stores the image ids.

- cellType:

  The column in colData that stores the cell types.

- ...:

  Any arguments passed into
  [`Kontextual`](https://sydneybiox.github.io/Statial/reference/Kontextual.md).

## Value

A data frame of original L values and Kontextual values evaluated over a
range of radii.

## Examples

``` r

data("kerenSCE")

kerenImage6 <- kerenSCE[, kerenSCE$imageID == "6"]

rsDf <- kontextCurve(
  cells = kerenSCE,
  from = "CD4_Cell",
  to = "Keratin_Tumour",
  parent = c("CD4_Cell", "Macrophages"),
  rs = seq(10, 510, 100),
  cores = 2
)
```
