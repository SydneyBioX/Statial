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
#' @param image A vector of images to subset the results to. If NULL we default to all images.
#' @param inhom  A logical value indicating whether to account for inhomogeneity.
#' @param edgeCorrect A logical value indicating whether to perform edge correction.
#' @param window Type of window for data, either `square`, `convex` or `concave`,
#'  passed into \code{\link[Statial]{makeWindow}}
#' @param window.length A tuning parameter for controlling the level of concavity when estimating concave windows.
#' Passed into \code{\link[Statial]{makeWindow}}
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
#' @importFrom data.table as.data.table 

Kontextual <- function(cells,
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
      spatialCoords = spatialCoords,
      image = image
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
      spatialCoords = spatialCoords,
      image = image
    )

    # cells <- mutate(cells, imageID = as.character(imageID))
    cells <- split(cells, cells$imageID)
  }

  if (!is(cells, "list")) {
    stop("Cells must be one of the following: SingleCellExperiment, SpatialExperiment, or a list of data.frames with imageID, cellType, and x and y columns")
  }
  
  
  
  # Specify cores for parrellel computing
  x <- runif(1) # nolint
  BPPARAM <- Statial:::.generateBPParam(cores = cores)
  
  r = sort(r)

  cells <- lapply(cells, function(image) {
    image$cellID <- factor(seq_len(nrow(image)))
    return(image)
  })

  images <- cells
  
  # Convert images to PPP objects
  imagesPPP <- bplapply(images, function(image){
    # make window
    ow <- Statial::makeWindow(image, window, window.length)
    
    # Define marks
    marks = data.frame(
      cellType = image$cellType,
      cellID = image$cellID
    )
    
    # Make PPP object
    imagePPP <- spatstat.geom::ppp(
      x = image$x,
      y = image$y,
      window = ow,
      marks = marks
    )
    
    return(imagePPP)
  }, BPPARAM = BPPARAM ) 
  

  # Calculate Areas
  Areas = bplapply(imagesPPP, function(image) {
    return(spatstat.geom::area(image))
  }, BPPARAM = BPPARAM )


  
  imagesInfo <- data.frame(
    imageID = names(images),
    images = I(images),
    imagesPPP = I(imagesPPP),
    area = I(Areas)
  )

  imagesInfoR <- expand_grid(
    imagesInfo,
    r = r) |> 
    mutate(id = paste(imageID, r, sep = "."))
  
  # Calculate close pairs for all images and radii
  closePairList <- bplapply(imagesPPP, function(imagePPP) {

    # Calculate close pairs
    closePairs <- spatstat.geom::closepairs(
      imagePPP, max(r, na.rm = TRUE),
      what = "ijd", distinct = FALSE
    ) |>
      data.frame()

    cellTypes <- imagePPP$marks$cellType
    names(cellTypes) <- imagePPP$marks$cellID
    closePairs$cellTypeI <- cellTypes[(closePairs$i)]
    closePairs$cellTypeJ <- cellTypes[(closePairs$j)]
    closePairs$i <- factor(closePairs$i, levels = imagePPP$marks$cellID)
    
    closePairs <- as.data.table(closePairs)
   
    # Subsetting close pairs if there is more R values.
    if(length(r) > 1){
      closePairList = lapply(r[-length(r)], function(r){
        return(closePairs[d < r])
      })
      
      closePairList = append(closePairList, list(closePairs))
      names(closePairList) = r
    } else{
      closePairList = list(closePairs)
      names(closePairList) = r
    }
    
    return(closePairList)
    
  }, BPPARAM = BPPARAM)
  
  
  # Adding closePairs to the image data.frame
  closePairList = unlist(closePairList, recursive = FALSE)
  closePairsDfs = data.frame(id = names(closePairList), 
                             closePairs = I(closePairList))
  imagesInfoR = left_join(imagesInfoR,
                          closePairsDfs,
                          by = join_by(id))
  

  
  # Perform edge correction (if needed)
  imagesInfoR$closePairs <- bpmapply(function(imagePPP, closePairs, r){
    
    if(edgeCorrect){
      edge <- .borderEdge(imagePPP, r)
      edge <- as.data.frame(edge)
      edge$i <- factor(imagePPP$marks$cellID, levels = imagePPP$marks$cellID)
      edge$edge <- 1 / edge$edge
    } else {
      edge = data.frame(
        i = factor(imagePPP$marks$cellID, levels = imagePPP$marks$cellID),
        edge = 1
      )
    }

    closePairs <- left_join(closePairs, edge[, c("i", "edge")], by = "i")
  }, 
  imagePPP = imagesInfoR$imagesPPP,
  closePairs = imagesInfoR$closePairs,
  r = imagesInfoR$r,
  SIMPLIFY = FALSE,
  BPPARAM = BPPARAM)


  # Create all combinations of specified parameters
  allCombinations <- expand_grid(
    parentDf,
    inhomL = inhom
  )

  # Create data frame for mapply
  kontextualDf <- merge(imagesInfoR, allCombinations, all = TRUE) |>
    mutate("test" = paste(from, "__", to, sep = ""))

  if ("parent_name" %in% names(kontextualDf)) {
    kontextualDf <- mutate(
      kontextualDf,
      "test" = paste(test, "__", parent_name, sep = "")
    )
  }

  # Calculate Kontextual values
  lVals <- bpmapply(
    KontextualCore,
    images = kontextualDf$images,
    r = kontextualDf$r,
    closePairs = kontextualDf$closePairs,
    from = kontextualDf$from,
    to = kontextualDf$to,
    parent = kontextualDf$parent,
    inhom = kontextualDf$inhomL,
    area = kontextualDf$area,
    SIMPLIFY = FALSE,
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
        "r"
      )
  } else {
    lValsClean <- lValsClean |>
      select(
        "imageID",
        "test",
        "original",
        "kontextual",
        "r",
        "inhomL"
      )
  }

  return(lValsClean)
}




#' @noRd
#'
#' @importFrom data.table dcast
KontextualCore <- function(images,
                           r,
                           from,
                           to,
                           parent,
                           closePairs,
                           area,
                           inhom = FALSE,
                           edge = FALSE,
                           includeOriginal = TRUE,
                           returnWeight = FALSE,
                           ...) {
  
  # Returns NA if to and from cell types not in image
  if (!(c(to, from) %in% unique(images$cellType) |> all())) {
    condL <- data.frame(original = NA, kontextual = NA)
    rownames(condL) <- paste(from, "__", to)
    return(condL)
  }
  
  
  child1 = from
  child2 = to
  
  # Convert closepairs to counts of cells next to other cells
  counts <- closePairs[, .(n = sum(edge)), by = .(i, cellTypeI, cellTypeJ)]
  counts <- dcast(counts, i + cellTypeI ~ cellTypeJ, value.var = "n", fill = 0)

  
  ########
  # Calculate statistics
  ########
  
  Kontextual <- .Kontext(closePairs, counts, child1, child2, parent, r, returnWeight)
  
  if (inhom) {
    L <- .Linhomfunction(closePairs, counts, child1, child2, r, area)
  } else {
    L <- .Lfunction(closePairs, counts, child1, child2, r, area)
  }
  
  if (returnWeight) {
    return(Kontext)
  }
  
  
  # return data frame of original and kontextual values.
  condL <- data.frame(
    original = L,
    kontextual = Kontextual
  )

  rownames(condL) <- paste(from, "__", to)


  return(condL)
}


#' @noRd
#'
#' @importFrom dplyr rename
validateDf <- function(cells, cellType, imageID, spatialCoords, image = NULL) {
  result <- try(
    {
      cells <- cells[, c(cellType, imageID, spatialCoords)] |>
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

  if (!is.null(image)) cells <- cells[cells$imageID %in% image, ]

  return(cells)
}




#' Test whether an object is a kontextualResult
#'
#' @param kontextualResult a object to test
#'
#' @examples
#' data <- data.frame()
#' if (!isKontextual(data)) print("Not a kontextualResult")
#'
#' @export isKontextual
#' @rdname isKontextual
isKontextual <- function(kontextualResult) {
  colNames <- c(
    "imageID",
    "test",
    "kontextual",
    "r"
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
#' kontextMat <- prepMatrix(CD4_Kontextual)
#'
#' @export prepMatrix
#' @rdname prepMatrix
#' @importFrom dplyr mutate
prepMatrix <- function(result,
                       replaceVal = 0,
                       column = NULL,
                       test = NULL) {
  mat <- NULL

  if ("kontextual" %in% colnames(result)) {
    mat <- result |>
      # Implement support for multiple values in other columns.
      dplyr::select(imageID, test, kontextual) |>
      tidyr::pivot_wider(names_from = test, values_from = kontextual, values_fill = replaceVal) |>
      tibble::column_to_rownames("imageID")
  }

  if ("primaryCellType" %in% colnames(result)) {
    if (!is.null(column)) {
      result$type <- result[, column]
    } else {
      result$type <- result$coef
    }

    mat <- result |>
      as.data.frame() |>
      # Implement support for multiple values in other columns.
      dplyr::mutate(
        test = paste(primaryCellType, otherCellType, marker, sep = "__")
      ) |>
      dplyr::select(imageID, test, type) |>
      tidyr::pivot_wider(
        names_from = test, values_from = type, values_fill = replaceVal
      ) |>
      tibble::column_to_rownames("imageID")
  }



  mat
}
