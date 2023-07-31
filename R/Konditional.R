#' Evaluation of pairwise cell relationships, conditional on a 3rd population.
#' 
#' @description 
#' Kontextual identifies the relationship between two cell types which are 
#' conditional on the spatial behaviour of a 3rd cell population, for a 
#' particular radius (r).
#' 
#'
#' @param cells A SingleCellExperiment, SpatialExperiment or a list of
#' data.frames containing columns specifying the imageID, cellType, and x and y
#' spatial coordinates.
#' @param parentDf A data frame from \code{\link[Statial]{parentCombinations}}
#' @param r Radii to evaluated pairwise relationships between from and to cells.
#' @param from The first cell type to be evaluated in the pairwise relationship.
#' @param to The second cell type to be evaluated in the pairwise relationship.
#' @param parent The parent population of the from cell type (must include from cell type).
#' @param image A vector of images to filter results to.
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
#' values along with the kontextual values.
#' @param spatialCoords The columns which contain the x and y spatial coordinates.
#' @param cellType The column which contains the cell types.
#' @param imageID The column which contains image identifiers.
#' @param cores Number of cores for parallel processing.
#' @return A kontextualResult object
#'
#' @examples
#' # Load data
#' data("kerenSCE")
#' 
#'
#' CD4_Kontextual <- Kontextual(
#'   cells = kerenSCE,
#'   r = 50,
#'   from = "Macrophages",
#'   to = "Keratin_Tumour",
#'   parent = c("Macrophages", "CD4_Cell"),
#'   image = "6"
#' )
#'
#'
#' head(CD4_Kontextual)
#'
#' @export Kontextual
#' @rdname Kontextual
#' @importFrom BiocParallel bpmapply
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble remove_rownames
#' @importFrom methods is
#' @importFrom stats runif

Kontextual <- function(cells,
                        r,
                        parentDf = NULL,
                        from = NULL,
                        to = NULL,
                        parent = NULL,
                        image = NULL,
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
  
  cells$imageID <- colData(cells)[,imageID]
  cells$cellType <- colData(cells)[,cellType]
  cellType <- "cellType"
  imageID <- "imageID"
  if(!is.null(image))cells <- cells[,cells$imageID %in% image]


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
      SummarizedExperiment::colData() |>
      data.frame()
  }

  if (is(cells, "SpatialExperiment")) {
    cells <- cbind(colData(cells), SpatialExperiment::spatialCoords(cells)) |>
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
  kontextualDf <- merge(imagesInfo, allCombinations, all = TRUE) |>
    mutate("test" = paste(from, "__", to, sep = ""))

  if ("parent_name" %in% names(kontextualDf)) {
    kontextualDf <- mutate(kontextualDf, "test" = paste(test, "__", parent_name, sep = ""))
  }

  x <- runif(1) # nolint

  BPPARAM <- .generateBPParam(cores = cores)

  # Calculate conditional L values
  lVals <- bpmapply(
    KontextualCore,
    images = kontextualDf$images,
    r = kontextualDf$r,
    from = kontextualDf$from,
    to = kontextualDf$to,
    parent = kontextualDf$parent,
    inhom = kontextualDf$inhom,
    edge = kontextualDf$edge,
    window = kontextualDf$window,
    window.length = kontextualDf$window.length,
    weightQuantile = kontextualDf$weightQuantile,
    includeZeroCells = kontextualDf$includeZeroCells,
    SIMPLIFY = FALSE,
    MoreArgs = list(includeOriginal = includeOriginal),
    BPPARAM = BPPARAM
  )

  # Combine data.frame rows
  lVals <- lVals |>
    bind_rows()

  lValsClean <- kontextualDf |>
    mutate("parent_name" = "") |>
    select(-c("images", "from", "to", "parent_name", "parent")) |>
    cbind(lVals) |>
    remove_rownames()

  if (includeOriginal == FALSE) {
    lValsClean <- lValsClean |>
      select(
        "imageID",
        "test",
        "kontextual",
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
        "kontextual",
        "r",
        "weightQuantile",
        "inhom",
        "edge",
        "includeZeroCells",
        "window",
        "window.length"
      )
  }
  
  return(lValsClean)
}




#' @noRd
#'
KontextualCore <- function(images,
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
  if (!(c(to, from) %in% unique(images$cellType) |> all())) {
    condL <- data.frame(original = NA, kontextual = NA)
    rownames(condL) <- paste(from, "__", to)
    return(condL)
  }

  # Calculated the child and parent values for the image
  kontextualVal <-
    inhomLParent(
      images,
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
    condL <- data.frame(kontextual = kontextualVal)
    return(condL)
  }

  originalVal <-
    inhomLParent(
      images,
      Rs = r,
      from = from,
      to = to,
      parent = unique(images$cellType),
      edgeCorrect = edge,
      inhom = inhom,
      weightQuantile = 1,
      original = TRUE,
      ...
    )

  # return data frame of original and kontextual values.
  condL <- data.frame(
    original = originalVal,
    kontextual = kontextualVal
  )

  return(condL)
}


#' @noRd
#'
#' @importFrom dplyr rename
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




#' Test whether an object is a kontextualResult
#'
#' @param kontextualResult a object to test
#'
#' @examples
#' data = data.frame()
#' if(!isKontextual(data)) print("Not a kontextualResult")
#'
#' @export isKontextual
#' @rdname isKontextual
isKontextual <- function(kontextualResult){
    
    colNames = c(
        "imageID",
        "test",
        "original",
        "kontextual",
        "r",
        "weightQuantile",
        "inhom",
        "edge",
        "includeZeroCells",
        "window",
        "window.length"
    )
    
    return(all(colNames %in% names(kontextualResult)))
}

#' Convert Kontextual or state changes result to a matrix for classification
#'
#' @param result a kontextual or state changes result data.frame.
#' @param replaceVal value which NAs are replaced with.
#' @param column The column which contains the scores that you want to select.
#' @param test A column containing which will be the column names of the expanded matrix.
#'
#' @examples
#' data("kerenSCE")
#' 
#'
#' CD4_Kontextual <- Kontextual(
#'   cells = kerenSCE,
#'   r = 50,
#'   from = "Macrophages",
#'   to = "Keratin_Tumour",
#'   parent = c("Macrophages", "CD4_Cell"),
#'   image = "6"
#' )
#'
#'
#' kontextMat = prepMatrix(CD4_Kontextual)
#'
#' @export prepMatrix
#' @rdname prepMatrix
#' @importFrom dplyr mutate
prepMatrix = function(result,
                      replaceVal = 0,
                      column = NULL,
                      test = NULL) {
  
  mat = NULL

  if("kontextual"%in%colnames(result)){
    
    mat <- result |> 
        # Implement support for multiple values in other columns.
        dplyr::select(imageID, test, kontextual) |> 
        tidyr::pivot_wider(names_from = test, values_from = kontextual, values_fill = replaceVal) |> 
        tibble::column_to_rownames("imageID") 
    
  }
  
  if("primaryCellType"%in%colnames(result)){
    
    if(!is.null(column)){
      result$type <- result[, column]
    }else{
      result$type <- result$coef
      }
    
    mat <- result |> 
      as.data.frame() |>
      # Implement support for multiple values in other columns.
      dplyr::mutate(test = paste(primaryCellType, otherCellType, marker, sep = "__")) |>
      dplyr::select(imageID, test, type) |> 
      tidyr::pivot_wider(names_from = test, values_from = type, values_fill = replaceVal) |> 
      tibble::column_to_rownames("imageID")

  }
    
    
  
  mat
}
    