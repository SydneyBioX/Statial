#' Test Kontextual relationships between conditions or with survival
#'
#' @description
#' For every triple `from` → `to` within `parent` (from [parentCombinations()]),
#' `kontextualTest()` tests whether the Kontextual relationship differs between
#' conditions, or is associated with survival, with patients as the units.
#'
#' The statistic is that of [Kontextual()]: each `from`-`to` pair within `r`
#' is weighted by the density of the parent population at the `from` cell over
#' that at the `to` cell (edge corrected), and averaged over the `from` cells in
#' proportion to the parent density around them. Instead of comparing these
#' values with zero, each image's value is compared with its exact expectation
#' and variance if the `to` cells were a random choice among the parent's cells.
#' A change in tissue composition between conditions therefore does not by
#' itself give a difference. The images' excesses over this expectation are
#' combined within patients and between conditions by a frailty GEE with a CR2
#' cluster-robust variance (Satterthwaite degrees of freedom), allowing for
#' clustering of the `to` labels within the parent, and by default adjusted for
#' the abundance of `from` (its log share of each image's cells). With a
#' survival outcome the test is a score test against the hazard.
#'
#' The excess is on the scale of Kontextual's K function (an area, in the
#' squared units of the coordinates): an image's [Kontextual()] value is
#' \eqn{\sqrt{K/\pi} - r}, and \eqn{K = \pi r^2} is expected when the `to` cells
#' are a random choice among the parent's cells.
#'
#' The computation is done by spicyR; the result is a `SpicyResults` object, so
#' [spicyR::topPairs()], [spicyR::signifPlot()] and [spicyR::spicyBoxPlot()]
#' work on it.
#'
#' @param cells A SingleCellExperiment, SpatialExperiment, data frame, or list
#'   of data frames (one per image) with a row per cell.
#' @param parentDf A data frame from [parentCombinations()].
#' @param condition The column with each image's condition (constant within a
#'   patient), or a `survival::Surv` column of survival outcomes.
#' @param r The radius.
#' @param subject The column identifying patients, when patients have several
#'   images.
#' @param covariates Image- or patient-level covariates to adjust for.
#' @param from,to Optional: keep only the triples with these cell types.
#' @param cellType,imageID,spatialCoords The columns of the cell type, image and
#'   coordinates.
#' @param window The window used for the edge correction: `"convex"` (the
#'   default, as in [Kontextual()]) or `"square"`.
#' @param edgeCorrect Whether to correct the parent densities at the edges of
#'   the window.
#' @param adjustAbundance Whether to adjust for the abundance of `from`.
#' @param labelClustering Whether to allow for clustering of the `to` labels
#'   within the parent.
#' @param ref The reference level of `condition`.
#'
#' @return A `SpicyResults` object with one row per triple (`from`, `to`,
#'   `parent`).
#'
#' @seealso [Kontextual()] for the per-image values.
#'
#' @examples
#' data("kerenSCE")
#'
#' kerenSCE$event <- 1 - kerenSCE$Censored
#' kerenSCE$survival <- survival::Surv(kerenSCE$Survival_days_capped, kerenSCE$event)
#'
#' parentDf <- parentCombinations(
#'   all = unique(kerenSCE$cellType),
#'   parentList = list(tcells = c("CD4_Cell", "CD8_Cell", "Tregs"))
#' )
#'
#' res <- kontextualTest(kerenSCE, parentDf, condition = "survival", r = 100,
#'                       from = "Keratin_Tumour")
#' spicyR::topPairs(res)
#'
#' @export
kontextualTest <- function(cells,
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
                           labelClustering = TRUE,
                           ref = NULL) {
  window <- match.arg(window)
  if (is.list(cells) && !is.data.frame(cells) && !methods::is(cells, "SummarizedExperiment")) {
    cells <- do.call(rbind, cells)
  }
  spicyR::kontextualEngine(
    cells = cells, parentDf = parentDf, condition = condition, subject = subject,
    covariates = covariates, r = r, from = from, to = to, imageID = imageID,
    cellType = cellType, spatialCoords = spatialCoords,
    adjustAbundance = adjustAbundance, labelClustering = labelClustering,
    edgeCorrect = edgeCorrect, window = if (window == "square") "rectangle" else "convex",
    ref = ref
  )
}
