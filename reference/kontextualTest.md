# Test Kontextual relationships between conditions or with survival

For every triple \`from\` → \`to\` within \`parent\` (from
\[parentCombinations()\]), \`kontextualTest()\` tests whether the
Kontextual relationship differs between conditions, or is associated
with survival, with patients as the units.

The statistic is that of \[Kontextual()\]: each \`from\`-\`to\` pair
within \`r\` is weighted by the density of the parent population at the
\`from\` cell over that at the \`to\` cell (edge corrected), and
averaged over the \`from\` cells in proportion to the parent density
around them. Instead of comparing these values with zero, each image's
value is compared with its exact expectation and variance if the \`to\`
cells were a random choice among the parent's cells. A change in tissue
composition between conditions therefore does not by itself give a
difference. The images' excesses over this expectation are combined
within patients and between conditions by a frailty GEE with a CR2
cluster-robust variance (Satterthwaite degrees of freedom), allowing for
clustering of the \`to\` labels within the parent, and by default
adjusted for the abundance of \`from\` (its log share of each image's
cells). With a survival outcome the test is a score test against the
hazard.

The excess is on the scale of Kontextual's K function (an area, in the
squared units of the coordinates): an image's \[Kontextual()\] value is
\\\sqrt{K/\pi} - r\\, and \\K = \pi r^2\\ is expected when the \`to\`
cells are a random choice among the parent's cells.

The computation is done by spicyR; the result is a \`SpicyResults\`
object, so \[spicyR::topPairs()\], \[spicyR::signifPlot()\] and
\[spicyR::spicyBoxPlot()\] work on it.

## Usage

``` r
kontextualTest(
  cells,
  parentDf,
  condition,
  r,
  subject = NULL,
  covariates = NULL,
  from = NULL,
  to = NULL,
  cellType = "cellType",
  imageID = "imageID",
  spatialCoords = c("x", "y"),
  window = c("convex", "square"),
  edgeCorrect = TRUE,
  adjustAbundance = TRUE,
  labelClustering = FALSE,
  ref = NULL
)
```

## Arguments

- cells:

  A SingleCellExperiment, SpatialExperiment, data frame, or list of data
  frames (one per image) with a row per cell.

- parentDf:

  A data frame from \[parentCombinations()\].

- condition:

  The column with each image's condition (constant within a patient), or
  a \`survival::Surv\` column of survival outcomes.

- r:

  The radius.

- subject:

  The column identifying patients, when patients have several images.

- covariates:

  Image- or patient-level covariates to adjust for.

- from, to:

  Optional: keep only the triples with these cell types.

- cellType, imageID, spatialCoords:

  The columns of the cell type, image and coordinates.

- window:

  The window used for the edge correction: \`"convex"\` (the default, as
  in \[Kontextual()\]) or \`"square"\`.

- edgeCorrect:

  Whether to correct the parent densities at the edges of the window.

- adjustAbundance:

  Whether to adjust for the abundance of \`from\`.

- labelClustering:

  Whether to allow for clustering of the \`to\` labels within the
  parent.

- ref:

  The reference level of \`condition\`.

## Value

A \`SpicyResults\` object with one row per triple (\`from\`, \`to\`,
\`parent\`).

## See also

\[Kontextual()\] for the per-image values.

## Examples

``` r
data("kerenSCE")

kerenSCE$event <- 1 - kerenSCE$Censored
kerenSCE$survival <- survival::Surv(kerenSCE$Survival_days_capped, kerenSCE$event)

parentDf <- parentCombinations(
  all = unique(kerenSCE$cellType),
  parentList = list(tcells = c("CD4_Cell", "CD8_Cell", "Tregs"))
)

res <- kontextualTest(kerenSCE, parentDf, condition = "survival", r = 100,
                      from = "Keratin_Tumour")
spicyR::topPairs(res)
#>                                            from       to parent hazard_ratio_sd
#> Keratin_Tumour__CD8_Cell__tcells Keratin_Tumour CD8_Cell tcells       0.7330915
#> Keratin_Tumour__CD4_Cell__tcells Keratin_Tumour CD4_Cell tcells       1.2586190
#>                                    p_value     p_adj
#> Keratin_Tumour__CD8_Cell__tcells 0.6532635 0.6895178
#> Keratin_Tumour__CD4_Cell__tcells 0.6895178 0.6895178
```
