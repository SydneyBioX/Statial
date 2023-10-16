
#' Cell permutation for Kontextual
#' 
#' @description 
#' Function which randomises specified cells in an image and calculates
#' the `Kontextual` value. This can be used to estimate the null distribution, 
#' of the parent cell population for significance testing.
#'
#' @param cells A single image data frame from a SingleCellExperiment object
#' @param nSim Number of randomisations which will be calculated.
#' @param r Radius to evaluated pairwise relationships between from and to cells.
#' @param from The first cell type to be evaluated in the pairwise relationship.
#' @param to The second cell type to be evaluated in the pairwise relationship.
#' @param parent The parent population of the from cell type (must include from cell type).
#' @param image A vector of images to subset the results to. If NULL we default to all images.
#' @param returnImages A logical value to indicate whether the function should
#' return the randomised images along with the Kontextual values.
#' @param inhom A logical value indicating whether to account for inhomogeneity.
#' @param edge A logical value indicating whether to perform edge correction.
#' @param cores Number of cores for parallel processing.
#' @param spatialCoords A character vector containing the names of the two
#'     spatial dimansions in the data. Defaults to `c("x", "y")`.
#' @param cellType The name of the cell type field in the data. Defualts to
#'     "cellType".
#' @param imageID The name of the image ID field in the data. Defualts to
#'     "imageID".
#' @param ... Any arguments passed into \code{\link[Statial]{Kontextual}}
#' @return A data frame containing Kontextual value for each randomised image.
#' If `returnImages = TRUE` function will return a list with Kontextual values
#' and the randomised images.
#'
#' @examples
#' data("kerenSCE")
#' 
#' kerenImage6 = kerenSCE[, kerenSCE$imageID =="6"]
#'
#' relabelResult <- relabelKontextual(
#'   cells = kerenImage6,
#'   nSim = 5,
#'   r = 250,
#'   from = "CD4_Cell",
#'   to = "Keratin_Tumour",
#'   parent = c("CD4_Cell", "Macrophages"),
#'   cores = 2
#' )
#'
#' @export
#' @rdname relabelKontextual
#' @importFrom dplyr mutate
#' @importFrom BiocParallel  bplapply
relabelKontextual <- function(cells,
                              nSim = 1,
                              r,
                              from,
                              to,
                              parent,
                              image = NULL,
                              returnImages = FALSE,
                              inhom = TRUE,
                              edge = FALSE,
                              cores = 1,
                              spatialCoords = c("x", "y"),
                              cellType = "cellType",
                              imageID = "imageID",
                              ...) {
  
  
  if (is(cells, "SingleCellExperiment")) {
    cells <- cells |>
      SingleCellExperiment::colData() |>
      data.frame()
  }
  
  if (is(cells, "SpatialExperiment")) {
    cells <- cbind(colData(cells), spatialCoords(cells)) |>
      data.frame()
  }
  
  if (is(cells, "data.frame")) {
    cells <- validateDf(
      cells,
      imageID = imageID,
      cellType = cellType,
      spatialCoords = spatialCoords,
      image = image
    )
  }
  
  if (!is(cells, "data.frame")) {
    stop("Cells must be one of the following: SingleCellExperiment, SpatialExperiment, or a list of data.frames with imageID, cellType, and x and y columns")
  }
  
  imageArray <- replicate(nSim, cells, simplify = FALSE)
  
  
  # relabel cells in parent population for all images in imageArray
  relabeled <- bplapply(imageArray,
                        relabel,
                        labels = parent,
                        BPPARAM = MulticoreParam(workers = cores)
  )
  
  relabeled <- c(list(cells), relabeled)
  names(relabeled) <- as.character(seq_along(relabeled))
  
  parentDf <- data.frame(
    from = from,
    to = to,
    parent = I(list(parent))
  )
  
  # Calculated the child and parent values for the relabeled images
  relabeledDf <- Kontextual(
    cells = relabeled,
    r = r,
    parentDf = parentDf,
    inhom = inhom,
    edgeCorrect = edge,
    cores = cores,
    ...
  )
  
  relabeledDf <- relabeledDf |>
    select("imageID", "original", "kontextual", "r") |>
    mutate(type = ifelse(imageID == 1, "original", "randomised"))
  
  if (returnImages) {
    return(list(relabeledDf, relabeled))
  }
  
  return(relabeledDf)
}



#' Permute all specified cells labels in a single image
#' 
#' @description 
#' This function relabels all specified cells within a single image, to 
#' estimate the null distribution of cell population specified.
#'
#' @param image A single image from a Single Cell Experiment object.
#' @param labels A vector of CellTypes labels to be permuted If NULL all cells
#' labels will be radomised.
#'
#' @return A data frame containing all pairwise cell relationships and their
#' corresponding parent
#'
#' @examples
#' data("kerenSCE")
#' 
#' kerenImage6 = kerenSCE[, kerenSCE$imageID =="6"]
#' 
#' kerenImage6 <- kerenImage6 |>
#'          SingleCellExperiment::colData() |>
#'          data.frame()
#'
#' # Permute CD8 T cells and T cell labels in the image
#' relabeledImage <- relabel(kerenImage6, labels = c("p53", "Keratin+Tumour"))
#' plot(relabeledImage)
#'
#' @export
#' @rdname relabelKontextual
#' @import dplyr
#'
relabel <- function(image, labels = NULL) {
  # if labels are NULL relabel the whole image, otherwise relabel just the specified marks
  if (is.null(labels)) {
    relabeledCells <- image |> mutate(cellType = sample(cellType))
  } else {
    # split up cells into subset which will be relabeled, and subset which wont be relabeled
    notRelabel <- image |> filter(!(cellType %in% labels))
    toRelabel <- image |> filter(cellType %in% labels)

    # relabel toRelabel cells
    relabeled <- toRelabel |> mutate(cellType = sample(cellType))

    # combine relabeled cells and non relabeled cells
    relabeledCells <- rbind(relabeled, notRelabel)
  }

  return(relabeledCells)
}
