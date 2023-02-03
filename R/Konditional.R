#' Evaluation of pairwise cell relationships, conditional on a 3rd population.
#' 
#' @description 
#' Konditional identifies the relationship between two cell types which are 
#' conditional on the spatial behaviour of a 3rd cell population, for a 
#' particular radius (r).
#' 
#'
#' @param cells A SingleCellExperiment, SpatialExperiment or a list of
#' data.frames containing columns specifying the imageID, cellType, and x and y
#' spatial coordinates.
#' @param parentDf A data frame from \code{\link[Statial]{parentCombinations}}
#' @param r Radius to evaluated pairwise relationships between from and to cells.
#' @param from The first cell type to be evaluated in the pairwise relationship.
#' @param to The second cell type to be evaluated in the pairwise relationship.
#' @param parent The parent population of the from cell type (must include from cell type).
#' @param inhom  A logical value indicating whether to account for inhomogeneity.
#' @param edgeCorrect A logical value indicating whether to perform edge correction.
#' @param window Type of window for data, either `square`, `convex` or `concave`,
#'  passed into \code{\link[Statial]{makeWindow}}
#' @param window.length A tuning parameter for controlling the level of concavity when estimating concave windows.
#' Passed into \code{\link[Statial]{makeWindow}}
#' @param weightQuantile A decimal value indicating what quantile of parent
#' density used to weight the `from` cells.
#' @param includeZeroCells A logical value indicating whether to include cells
#' with zero counts in the pairwise association calculation.
#' @param includeOriginal A logical value to return the original L function
#' values along with the konditional values.
#' @param spatialCoords The columns which contain the x and y spatial coordinates.
#' @param cellType The column which contains the cell types.
#' @param imageID The column which contains image identifiers.
#' @param cores Number of cores for parallel processing.
#' @return A konditionalResult object
#'
#' @examples
#' # Load data
#' data("headSCE")
#'
#' CD4_Konditional <- Konditional(
#'   cells = headSCE,
#'   r = 50,
#'   from = "TC_CD4",
#'   to = "SC5",
#'   parent = c("TC_CD4", "TC_CD8"),
#'   cores = 40
#' )
#'
#'
#' head(CD4_Konditional)
#'
#' @export Konditional
#' @rdname Konditional
#' @import dplyr
#' @import tidyr
#' @import BiocParallel
#' @import SingleCellExperiment
#' @importFrom tibble remove_rownames
#' @importFrom methods is
#' @importFrom stats runif

Konditional <- function(cells,
                        parentDf = NULL,
                        r,
                        from = NULL,
                        to = NULL,
                        parent = NULL,
                        inhom = TRUE,
                        edgeCorrect = FALSE,
                        window = "convex",
                        window.length = NA,
                        weightQuantile = .80,
                        includeZeroCells = TRUE,
                        includeOriginal = TRUE,
                        spatialCoords = c("x", "y"),
                        cellType = "cellType",
                        imageID = "imageID",
                        cores = 1) {
  if (is.null(parentDf) &
    is.null(from) &
    is.null(to) &
    is.null(parent)) {
    stop("Please specificy a parentDf (obtained from parentCombinations), or from, to, and parent cellTypes")
  } else if (!is.null(from) &
    !is.null(to) &
    !is.null(parent)
  ) {
    parentDf <- data.frame(
      from = from,
      to = to,
      parent = I(list(parent))
    )
  }


  if (is(cells, "list")) {
    cells <- lapply(
      cells,
      validateDf,
      imageID = imageID,
      cellType = cellType,
      spatialCoords = spatialCoords
    )
  }

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
      spatialCoords = spatialCoords
    )

    cells <- mutate(cells, imageID = as.character(imageID))
    cells <- split(cells, cells$imageID)
  }

  if (!is(cells, "list")) {
    stop("Cells must be one of the following: SingleCellExperiment, SpatialExperiment, or a list of data.frames with imageID, cellType, and x and y columns")
  }



  images <- cells

  imagesInfo <- data.frame(
    imageID = names(images),
    images = I(images)
  )

  # Create all combinations of specified parameters
  allCombinations <- expand_grid(
    parentDf,
    r = r,
    inhom = inhom,
    edge = edgeCorrect,
    window = window,
    window.length = window.length,
    weightQuantile = weightQuantile,
    includeZeroCells = includeZeroCells
  )

  # Create data frame for mapply
  konditionalDf <- merge(imagesInfo, allCombinations, all = TRUE) |>
    mutate("test" = paste(from, "__", to, sep = ""))

  if ("parent_name" %in% names(konditionalDf)) {
    konditionalDf <- mutate(konditionalDf, "test" = paste(test, "__", parent_name, sep = ""))
  }

  x <- runif(1) # nolint

  BPPARAM <- .generateBPParam(cores = cores)

  # Calculate conditional L values
  lVals <- bpmapply(
    KonditionalCore,
    image = konditionalDf$images,
    r = konditionalDf$r,
    from = konditionalDf$from,
    to = konditionalDf$to,
    parent = konditionalDf$parent,
    inhom = konditionalDf$inhom,
    edge = konditionalDf$edge,
    window = konditionalDf$window,
    window.length = konditionalDf$window.length,
    weightQuantile = konditionalDf$weightQuantile,
    includeZeroCells = konditionalDf$includeZeroCells,
    SIMPLIFY = FALSE,
    MoreArgs = list(includeOriginal = includeOriginal),
    BPPARAM = BPPARAM
  )

  # Combine data.frame rows
  lVals <- lVals |>
    bind_rows()

  lValsClean <- konditionalDf |>
    mutate("parent_name" = "") |>
    select(-c("images", "from", "to", "parent_name", "parent")) |>
    cbind(lVals) |>
    remove_rownames()

  if (includeOriginal == FALSE) {
    lValsClean <- lValsClean |>
      select(
        "imageID",
        "test",
        "konditional",
        "r",
        "weightQuantile",
        "inhom",
        "edge",
        "includeZeroCells",
        "window",
        "window.length"
      )
  } else {
    lValsClean <- lValsClean |>
      select(
        "imageID",
        "test",
        "original",
        "konditional",
        "r",
        "weightQuantile",
        "inhom",
        "edge",
        "includeZeroCells",
        "window",
        "window.length"
      )
  }
  
  #lValsClean <- methods::new("konditionalResult", lValsClean)
  
  return(lValsClean)
}




#' @noRd
#'
#' @import tidyverse
KonditionalCore <- function(image,
                            r,
                            from,
                            to,
                            parent,
                            inhom = TRUE,
                            edge = FALSE,
                            includeOriginal = TRUE,
                            weightQuantile = .80,
                            ...) {
  # Returns NA if to and from cell types not in image
  if (!(c(to, from) %in% unique(image$cellType) |> all())) {
    condL <- data.frame(original = NA, konditional = NA)
    rownames(condL) <- paste(from, "__", to)
    return(condL)
  }

  # Calculated the child and parent values for the image
  konditionalVal <-
    inhomLParent(
      image,
      Rs = r,
      from = from,
      to = to,
      parent = parent,
      edgeCorrect = edge,
      inhom = inhom,
      weightQuantile = weightQuantile,
      ...
    )


  if (includeOriginal == FALSE) {
    condL <- data.frame(konditional = konditionalVal)
    return(condL)
  }

  originalVal <-
    inhomLParent(
      image,
      Rs = r,
      from = from,
      to = to,
      parent = unique(image$cellType),
      edgeCorrect = edge,
      inhom = inhom,
      weightQuantile = 1,
      ...
    )

  # return data frame of original and konditional values.
  condL <- data.frame(
    original = originalVal,
    konditional = konditionalVal
  )

  return(condL)
}


#' @noRd
#'
#' @import tidyverse
validateDf <- function(cells, cellType, imageID, spatialCoords) {
  if (!("imageID" %in% names(cells)) ||
    !("cellType" %in% names(cells)) ||
    !("x" %in% names(cells)) ||
    !("y" %in% names(cells))) {
    result <- try(
      {
        cells <- cells |>
          rename(
            "cellType" = cellType,
            "imageID" = imageID,
            "x" = spatialCoords[1],
            "y" = spatialCoords[2]
          )
      },
      silent = TRUE
    )

    if (is(result, "try-error")) {
      stop("Please specifiy imageID or cellType or spatialCoords")
    }
  }
  return(cells)
}




#' @noRd
#'
#' @import tidyverse
isKonditional <- function(konditionalResult){
    
    colNames = c(
        "imageID",
        "test",
        "original",
        "konditional",
        "r",
        "weightQuantile",
        "inhom",
        "edge",
        "includeZeroCells",
        "window",
        "window.length"
    )
    
    return(all(colNames %in% names(konditionalResult)))
}
