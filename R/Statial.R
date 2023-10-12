#' @noRd
#'
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
preProcessing <- function(SCE, intensities) {
  
  
  #Check if cellID exists in image already - if duplicated, reassign to
  #ImageCellID
  if(any(c("cellID", "CellID") %in% colnames(colData(SCE)))) {
    if(any(duplicated(colData(SCE)$cellID))) {
      colnames(colData(SCE))[which(names(colData(SCE)) == "cellID")] <- "imageCellID"
    } 
    #TODO:Need to change CellID to cellID
    if(any(duplicated(colData(SCE)$CellID))) {
      colnames(colData(SCE))[which(names(colData(SCE)) == "CellID")] <- "ImageCellID"      
    } else {
      colnames(colData(SCE))[which(names(colData(SCE)) == "CellID")] <- "cellID"      
    }
  }
  
  if(!c("cellID") %in% colnames(colData(SCE))) {
    colData(SCE)$cellID <- rownames(colData(SCE))
  }
  
  
  intensitiesData <- data.frame(t(assay(SCE, intensities)))
  # spatialData <- data.frame(colData(SCE))
  
  # singleCellData <- cbind(spatialData[rownames(intensitiesData), ], intensitiesData)
  intensitiesData <- intensitiesData %>%
    mutate_all(~ ifelse(is.na(.), 0, .)) %>% #replace NAs with 0
    mutate_if(is.factor, as.character)
  intensitiesData <- intensitiesData[match(rownames(colData(SCE)), rownames(intensitiesData)), ]
  if (any(is.na(intensitiesData))) {
    stop("Number of cells in Intensities assay does not match number of cells in SCE")
  }
  colData(SCE) <- cbind(colData(SCE), intensitiesData) 
  
  # Identify factor columns
  
  factor_cols <- which(sapply(colData(SCE), is.factor))
  
  
  for(x in factor_cols) {
    colData(SCE)[, x] <- colData(SCE)[, x] |> as.character()
  }

  return(SCE)
}

#TODO: unify things that you don't have to type imageID = "slide" every single
#time. 

#' Calculate pairwise distance between cell types
#'
#' Calculates the euclidean distance from each cell to the nearest cell of each
#' type for a single image
#'
#' @param data the single cell data of interest
#' @param maxDist Maximum distance between pairs of points to be counted as close
#'   pairs.
#' @param distFun How to merge duplicate entries.
#'
#' @rdname distanceCalculator
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr
#'   select mutate inner_join arrange group_by slice rename relocate full_join
#' @importFrom spatstat.geom owin ppp closepairs
#' @importFrom tidyr pivot_wider
distanceCalculator <- function(data, maxDist = 200, distFun = "min") {
 
  ow <- spatstat.geom::owin(
    xrange = range(data$x),
    yrange = range(data$y)
  )
  pppData <- spatstat.geom::ppp(
    x = data$x,
    y = data$y,
    window = ow,
    marks = data$cellType
  )
  
  closePairData <- spatstat.geom::closepairs(pppData, rmax = maxDist, what = "ijd")
  distanceData <- data.frame(
    cellID = data$cellID[closePairData$i],
    cellType = data$cellType[closePairData$j],
    d = closePairData$d
  )
  
  if(distFun == "min") values_fn <- function(x){min(c(x,maxDist), na.rm = TRUE)}
  if(distFun == "abundance") values_fn <- function(x){sum(c(x,0)>0, na.rm = TRUE)}
  
  values_fill = values_fn(NULL)
  
  distanceData <- tidyr::pivot_wider(distanceData, names_from = cellType, values_from = d, 
                                     values_fn = values_fn, values_fill = values_fill)
  distanceData <- dplyr::left_join(data[,"cellID", drop = FALSE], distanceData, by = "cellID")|>
    tibble::column_to_rownames("cellID") 
  
  #distanceData[is.na(distanceData)] <- maxDist
  distanceData
  
  }


#' Wrapper to calculate pairwise distance between cell types by image
#'
#' Calculates the euclidean distance from each cell to the nearest cell of each
#' type
#'
#' @param cells
#'   A dataframe with a cellType column as well as x and y spatial coordinates.
#'   The dataframe must contain a imageID column and cellID (unique cell
#'   identifier's) column as well
#' @param maxDist
#'   The maximum distance considered.
#' @param imageID The name of the colData column that stores in the image ID.
#' @param spatialCoords The columns that store the spatial coordinates.
#' @param cellType The name of the colData column that stores the cell types.
#' @param redDimName The name of the reduced dimension to store the distances in.
#' @param distFun What distance function to use. Can be min or abundance.
#' @param nCores Number of cores for parallel processing.
#' @examples
#' data("kerenSCE")
#'
#' kerenSCE <- getDistances(kerenSCE,
#'   maxDist = 200
#' )
#'
#' @export
#' @rdname getDistances
#' @importFrom dplyr
#'   bind_rows rename_with contains starts_with vars across
#' @importFrom spatstat.geom owin ppp closepairs
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom purrr reduce
#' @importFrom stringr str_replace
#' @importFrom SummarizedExperiment colData assayNames
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr select contains mutate
#' @importFrom magrittr %>%
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
getDistances <- function(cells,
                         maxDist = NULL,
                         imageID = "imageID",
                         spatialCoords = c("x", "y"),
                         cellType = "cellType",
                         redDimName = "distances",
                         distFun = "min",
                         nCores = 1) {
  x <- runif(1)
  BPPARAM <- .generateBPParam(cores = nCores)
  #TODO: Make all functions accept dataframe + spatialExperiment with spatialCoords in spatialCoords slot
  if(!is(cells, "SingleCellExperiment"))stop("Currently this only accepts SpatialExperiment or SingleCellExperiment")
  
  cd <- as.data.frame(SingleCellExperiment::colData(cells))
  
  if(!any(c(cellType, imageID, spatialCoords)%in%colnames(cd))) stop("Either imageID, cellType or spatialCoords is not in your colData")
  
  cd <- cd[, c(cellType, imageID, spatialCoords)]
  colnames(cd) <- c("cellType", "imageID", "x", "y")
  if(is.null(colnames(cells))) colnames(cells) <- seq_len(ncol(cells))
  cd$cellID <- colnames(cells)
  

  cdFilt <- cd
  
  if(is(cdFilt$imageID, "factor"))cdFilt$imageID <- droplevels(cdFilt$imageID)
  
  
  if(is.null(maxDist))maxDist <- max(diff(range(cdFilt$x)), diff(range(cdFilt$y)))
  
  distances <- cdFilt |>
    split(~imageID) |>
    BiocParallel::bplapply(distanceCalculator,
                           maxDist = maxDist,
                           distFun = distFun,
                           BPPARAM = BPPARAM
    )
  
  
  
  distances <- distances |>
    dplyr::bind_rows() 
  
  SingleCellExperiment::reducedDim(cells, redDimName) <- distances[colnames(cells),]
  
  cells
}






#' Wrapper to calculate imhomogenous K function between a cell and surrounding
#' types on each image
#'
#' Calculate the imhomogenous K function (a measure of cell type abundance) for
#' each cell to other cell types
#'
#' @param cells
#'   A dataframe with a cellType column as well as x and y spatial coordinates.
#'   The dataframe must contain a imageID column and cellID (unique cell
#'   identifier's) column as well
#' @param r
#'   Radius to include in that calculation of pairwise abundance (K-function)
#'   between cells (can be a numeric or vector of radii)
#' @param distFun What distance function to use.
#' @param redDimName Name of the reduced dimension to store in sce.
#' @param cellType The name of the column in colData that stores the cell types.
#' @param imageID The name of the column in colData that Stores the image ids.
#' @param spatialCoords The names of the columns in colData that store the spatial coordinates.
#' @param nCores Number of cores for parallel processing
#' 
#'
#' @examples
#' library(dplyr)
#' data("kerenSCE")
#' 
#' singleCellDataCounts <- getAbundances(kerenSCE,
#'   r = 200,
#' )
#'
#' @export
#' @rdname getAbundances
#' @importFrom dplyr bind_rows
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData assayNames
#' @importFrom BiocParallel bplapply
getAbundances <- function(cells,
                          r = 200,
                          distFun = "abundance",
                          redDimName = "abundances",
                          cellType = "cellType",
                          imageID = "imageID",
                          spatialCoords = c("x","y"),
                          nCores = 1) {
  
  x <- runif(1)
  BPPARAM <- .generateBPParam(cores = nCores)
  
  if(!is(cells, "SingleCellExperiment"))stop("Currently this only accepts SpatialExperiment or SingleCellExperiment")
  
  cd <- as.data.frame(SingleCellExperiment::colData(cells))
  
  if(!any(c(cellType, imageID, spatialCoords)%in%colnames(cd))) stop("Either imageID, cellType or spatialCoords is not in your colData")
  
  cd <- cd[, c(cellType, imageID, spatialCoords)]
  colnames(cd) <- c("cellType", "imageID", "x", "y")
  if(is.null(colnames(cells))) colnames(cells) <- seq_len(ncol(cells))
  cd$cellID <- colnames(cells)
  
  cdFilt <- cd
  if(is(cdFilt$imageID, "factor"))cdFilt$imageID <- droplevels(cdFilt$imageID)
  maxDist <- r

  if(is.null(maxDist))maxDist <- max(diff(range(cdFilt$x)), diff(range(cdFilt$y)))
  
  distances <- cdFilt |>
    split(~imageID) |>
    BiocParallel::bplapply(distanceCalculator,
                           maxDist = maxDist,
                           distFun = distFun,
                           BPPARAM = BPPARAM
    )
  
  
  
  distances <- distances |>
    dplyr::bind_rows() 
  
  SingleCellExperiment::reducedDim(cells, redDimName) <- distances[colnames(cells),]
  
  cells
}


#' Calculate the level of marker contamination of each cell
#'
#' Calculates contamination scores using a random forest classification
#'
#' @param cells
#'   A SingleCellExperiment or SpatialExperiment with a cellType column as well as marker intensity information
#'   corresponding to each cell. 
#' @param markers
#'   A vector of markers that proxy a cell's state. If NULL, all markers 
#'   will be used.
#' @param num.trees Number of trees to be used in the random forest classifier
#' @param verbose
#'   A logical indicating whether information about the final random forest
#'   model should be outputted.
#' @param missingReplacement
#'   A default value to replace missing marker intensities for classification.
#' @param assay The assay in the sce that contains the marker expressions.
#' @param cellType The name of the column in colData that stores the cell types.
#' @param redDimName The redDimName to store the output in the sce.
#'
#' @examples
#' data("kerenSCE")
#'
#' singleCellDataDistancesContam <- calcContamination(
#'   kerenSCE
#' )
#'
#' @export
#' @rdname calcContamination
#' @importFrom dplyr mutate select bind_rows across all_of
#' @importFrom SummarizedExperiment colData assayNames
#' @importFrom ranger ranger
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom stringr str_replace
#' @importFrom SingleCellExperiment reducedDim reducedDim<-
calcContamination <- function(cells,
                              markers = NULL,
                              num.trees = 100,
                              verbose = FALSE,
                              missingReplacement = 0,
                              assay = "intensities",
                              cellType = "cellType",
                              redDimName = "contaminations"
) {
  
  singleCellData <- cells
  if(is.null(markers)) {
    markers <- rownames(singleCellData)
  }
  
  
  if(!all(markers %in% colnames(colData(singleCellData)))) {
    singleCellDataClean <- singleCellData |>
      preProcessing(intensities = assay)
  } else {
    singleCellDataClean <- singleCellData
  }
  singleCellData <- as.data.frame(colData(singleCellDataClean))  
  
  
  
  rfData <- singleCellData |>
    dplyr::select(cellType, all_of(markers)) |>
    dplyr::mutate(dplyr::across(any_of(markers), function(x) ifelse(is.nan(x) | is.na(x), 0, x)))
  
  rfData$cellType <- rfData[[cellType]]
  
  rfModel <- ranger::ranger(
    as.factor(cellType) ~ .,
    data = rfData,
    num.trees = num.trees,
    probability = TRUE
  )
  
  if (verbose == TRUE) {
    print(rfModel)
  }
  
  predictions <- predict(rfModel, rfData)$predictions
  
  maxn <- function(n) function(x) order(x, decreasing = TRUE)[!is.na(x)][n]
  
  rfData <- cbind(rfData, predictions) 
  
  #check for duplicates - this can happen if your cells are annotated "CD4" for 
  #"CD4 T cells" and you also have a gene marker "CD4"
  if(!is.empty(which(duplicated(names(rfData))))) {
    duplicate_columns <- which(duplicated(names(rfData)))
    print(paste("You have duplicate columns at columns ", paste0(duplicate_columns, collapse = ",")))
  }
  
  rfData <- rfData %>%
    dplyr::mutate(
      rfMaxCellProb = apply(
        .[colnames(predictions)],
        1,
        function(x) x[maxn(1)(x)]
      )
    ) 
  rfData <- rfData %>%
    dplyr::mutate(
      rfSecondLargestCellProb = apply(
        .[colnames(predictions)],
        1,
        function(x) x[maxn(2)(x)]
      )
    )
  rfData <- rfData %>%
    dplyr::mutate(
      rfMainCellProb = apply(
        .[c("cellType", colnames(predictions))],
        1,
        function(x) as.numeric(x[x["cellType"]])
      )
    ) %>%
    #dplyr::select(-colnames(predictions)) %>%
    tibble::rownames_to_column("cellID") %>%
    dplyr::mutate(cellID = stringr::str_replace(cellID, "cellID", "")) %>%
    select(-any_of(markers))
  
  # singleCellData2 <- singleCellData %>%
  #   dplyr::left_join(rfData, by = c("cellID"))

  redDim <- rfData
  rownames(redDim) <- colnames(cells)
  
  
  reducedDim(cells, redDimName) <- redDim
  
  return(cells)
}



#' First layer wrapper function to build linear models measuring state changes
#'
#' Builds linear models measuring marker based state changes in a cell type
#' based of the proximity or abundance of another cell type. The function
#' provides the option to build robust and mixed linear model variants
#'
#'
#'
#' @param cells
#'   A dataframe with a imageID, cellType, and marker intensity column along
#'   with covariates (e.g. distance or abundance of the nearest cell type) to
#'   model cell state changes
#' @param marker
#'  A vector of markers that proxy a cell's state. If NULL, all markers 
#'  will be used.
#' @param from A vector of cell types to use as the primary cells. If NULL,
#'  all cell types  will be used.
#' @param to A vector of cell types to use as the interacting cells. If NULL,
#'  all cell types  will be used.
#' @param image A vector of images to filter to. If null all images will be used.
#' @param type What type of state change. This value should be in reduced dimensions.
#' @param assay The assay in the sce that contains the marker expressions.
#' @param cellType The column in colData that stores the cell types.
#' @param imageID The column in colData that stores the image ids.
#' @param contamination If TRUE, use the contamination scores that have previously
#'  been calculate. Otherwise a name of which reduced dimension contains the scores.
#' @param minCells The minimum number of cells required to fit a model.
#' @param verbose A logical indicating if messages should be printed
#' @param timeout
#'   A maximum time allowed to build each model. Setting this may be important
#'   when building rlm mixed linear models
#' @param nCores Number of cores for parallel processing
#'
#' @examples
#' library(dplyr)
#' data("kerenSCE")
#' 
#' kerenSCE <- kerenSCE[, kerenSCE$imageID %in% c(5,6)]
#'
#' kerenSCE <- getDistances(kerenSCE,
#'   maxDist = 200,
#' )
#'
#' imageModels <- calcStateChanges(
#'   cells = kerenSCE,
#'   from = "Macrophages",
#'   to = "Tumour"
#' )
#'
#' @export
#' @rdname calcStateChanges
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom dplyr
#'   arrange group_by  summarise_at mutate bind_rows left_join filter left_join
#' @importFrom tidyr gather
#' @importFrom SummarizedExperiment colData assay
calcStateChanges <- function(cells,
                            marker = NULL,
                            from = NULL,
                            to = NULL,
                            image = NULL,
                            type = "distances",
                            assay = 1,
                            cellType = "cellType",
                            imageID = "imageID",
                            contamination = NULL,
                            minCells = 20,
                            verbose = FALSE,
                            timeout = 10,
                            nCores = 1) {
  
  if(is.null(marker)) {
    marker <- rownames(cells)
  }
  
  if(!is.null(image)){
    cells <- cells[,colData(cells)[,imageID]%in%image]
  }
  
  if(is.null(to)) {
    to <- unique(colData(cells)[,cellType])
    if(class(to) == "factor") {
      to <- droplevels(to)
    }
  }
  
  if(is.null(from)) {
    from <- to
  }
  
  if(!is.null(contamination)){
    if(contamination == TRUE) contamination = "contaminations"
  }
  
  cells <- cells[, colData(cells)[,cellType]%in%from]
  distances <- SingleCellExperiment::reducedDim(cells, type)
  distances <- distances[, to, drop = FALSE]
  intensities <- as.data.frame(t(SummarizedExperiment::assay(cells, assay)))
  intensities <- intensities[,marker,drop = FALSE]
  
  splitDist <- split(distances, ~ colData(cells)[, imageID] + colData(cells)[, cellType], sep = "51773")
  splitInt <- split(intensities, ~ colData(cells)[, imageID] + colData(cells)[, cellType], sep = "51773")
  if(is.null(contamination)) {
    contaminations <- data.frame(madeUp = rep(-99, ncol(cells)))
  }else{
    contaminations <- SingleCellExperiment::reducedDim(cells, contamination)
    contaminations <- dplyr::select(contaminations, -cellID, -cellType, -rfMaxCellProb, -rfSecondLargestCellProb, -rfMainCellProb)
  }
  
    splitCon <- split(contaminations, ~ colData(cells)[, imageID] + colData(cells)[, cellType], sep = "51773")
  
  use <- unlist(lapply(splitDist, nrow)) > minCells
  x <- runif(1)
  BPPARAM <- .generateBPParam(cores = nCores)
  
  
  allModels <- bpmapply(calculateChangesMarker, 
                      distances =  splitDist[use], 
                      intensities = splitInt[use],
                      contaminations = splitCon[use],
                      BPPARAM = BPPARAM,
                      SIMPLIFY = FALSE
  )
  
  allModels <- dplyr::bind_rows(allModels, .id = "tmp")
  nam <- strsplit(allModels$tmp, 51773)
  allModels$primaryCellType <- unlist(lapply(nam,function(x)x[2]))
  allModels$imageID <- unlist(lapply(nam,function(x)x[1]))
  allModels$fdr <- p.adjust(allModels$pval, "fdr")
  df <- dplyr::select(allModels, imageID, primaryCellType, otherCellType, marker, coef, tval, pval, fdr)
  df[order(df$pval),]
}




#' @importFrom limma lmFit
calculateChangesMarker <- function(distances, intensities, contaminations, nCores){
  test <- apply(distances, 2, function(x){
    if(length(unique(x))>1){
      if(contaminations[1,1] == -99){design <- data.frame(coef = 1, cellType = x)
      }else{
        contaminations <- contaminations[, !is.na(colSums(contaminations)),drop = FALSE]
        design <- data.frame(coef = 1, cellType = x, contaminations[,-ncol(contaminations)])
      }
      exprs <- t(intensities)
      fit <- .quiet(limma::lmFit(exprs, design, trend = rep(1,nrow(exprs)), verbose = FALSE))
      
      df <- data.frame(
      marker = colnames(intensities),
      coef = fit$coef[,"cellType"],
      tval = (fit$coef/fit$stdev.unscaled/fit$sigma)[,"cellType"]
      )
      df$pval <- 2 * pt(-abs(df$t), df = fit$df.residual)
      rownames(df) <- NULL
    
      df
    }
    }, simplify = FALSE)
    
  test <- dplyr::bind_rows(test, .id = "otherCellType")
  
  
}
                      
                        
.quiet <- function (x, print_cat = TRUE, message = TRUE, warning = TRUE) 
{
  stopifnot(is.logical(print_cat) && length(print_cat) == 1)
  stopifnot(is.logical(message) && length(message) == 1)
  stopifnot(is.logical(warning) && length(warning) == 1)
  if (print_cat) 
    sink(tempfile(), type = "out")
  on.exit(if (print_cat) sink())
  if (warning && message) 
    invisible(force(suppressMessages(suppressWarnings(x))))
  else if (warning && !message) 
    invisible(suppressWarnings(force(x)))
  else if (!warning && message) 
    invisible(suppressMessages(force(x)))
  else invisible(force(x))
}                       
                        
                        




### Section 1: LM Calculations ##############################################################




### Section 4:Visualise Image ##############################################################

#' Visualise Cell-Cell Marker Relationships
#'
#' Helper functions to visualise OLS model fits for image based state models
#'
#'image,
#'
#' @param cells
#'   A SingleCellExperiment that has had distances already calculated.
#' @param type The name of the reduced dimension to use for the x-axis.
#' @param image An image to subset to.
#' @param imageID
#'   Identifier name of the image in the imageID column to be visualised
#' @param from
#'   A character indicating the name of the cell type (from the cellType column) whose cell state is being investigated in
#' @param to
#' A character indicating the name of the cell type (from the cellType column) who may be influencing the cell state of another cell type
#' @param marker
#'   The marker of interest.
#' @param assay Name of the assay that stores the marker expression.
#' @param cellType The name of the column in colData that stores the cell types.
#' @param imageID The name of the column in colData that stores the image ids.
#' @param spatialCoords The names of the columns in colData that store the spatial coordinates.
#' @param size
#'   Aesthetic numerical variable determining the size of the displayed cells
#' @param shape
#'   Aesthetic variable determining the shape grouping of the displayed cells
#' @param interactive
#'   Logical indicating if the output visualisation should be a interactive (plotly)
#' @param plotModelFit
#'   Logical indicating if fitted values should be plotted or actual intensities for marker specified. The default is to plot actual intensities
#' @param method 
#'   The method to build the model with. Currently the only option is "lm". However, capabilities may be expanded in the future
#'
#' @examples
#' library(dplyr)
#' data("kerenSCE")
#'
#' kerenSCE <- getDistances(kerenSCE)
#' 
#' p <- plotStateChanges(
#'   cells = kerenSCE,
#'   type = "distances",
#'   image = "6",
#'   from = "Keratin_Tumour",
#'   to = "Macrophages",
#'   marker = "p53",
#'   size = 1,
#'   shape = 19,
#'   interactive = FALSE,
#'   plotModelFit = FALSE,
#'   method = "lm")
#' 
#' p
#'
#' @export
#' @rdname plotStateChanges
#' @importFrom dplyr filter left_join
#' @importFrom ggplot2
#'   ggplot scale_fill_distiller stat_density_2d geom_point theme_classic
#'   aes_string ggtitle facet_wrap aes xlab ylab ggtitle autoplot scale_colour_gradientn
#' @importFrom plotly ggplotly
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
#' @importFrom SingleCellExperiment reducedDimNames
plotStateChanges <- function(cells,
                                       image,
                                       from,
                                       to,
                                       marker,
                                       type = "distances",
                                       assay = 1,
                                       cellType = "cellType",
                                       imageID = "imageID",
                                       spatialCoords = c("x","y"),
                                       size = 1,
                                       shape = 19,
                                       interactive = FALSE,
                                       plotModelFit = FALSE,
                                       method = "lm") {
  

  if(!imageID %in% colnames(colData(cells))) {
    stop("The argument imageID needs to exist in the data")
  }

  if(!image %in% colData(cells)[,imageID]) {
    stop("The image needs to exist in the data")
  }
  
  cells <- cells[, colData(cells)[,imageID] %in% image]
  
  if(!type %in% reducedDimNames(cells)) {
    stop("The reduced dimension needs to exist in the data")
  }
  
  
  data <- data.frame(t(assay(cells, assay)), reducedDim(cells, type), colData(cells)) 
  
  data$imageID <- data[, imageID]
  data$cellType <- data[, cellType]
  
  # data <- metadata(SCE)[[metadata_name]]
  
  if(!marker %in% colnames(data)) {
    stop("The marker needs to exist in the data")
  }

  if(!cellType %in% colnames(data)) {
    stop("cellType needs to exist in the colData")
  }
  
  if(!all(spatialCoords %in% colnames(data))) {
    stop("spatialCoords needs to exist in the colData")
  }
  
  data$x <- data[, spatialCoords[1]]
  data$y <- data[, spatialCoords[2]]
  
  data$marker <- data[, marker, drop = TRUE]
  data$fittedValues <- NA
  
  contams <- reducedDim(cells, "contaminations")
  names(contams)[names(contams) == "rfMainCellProb"] <- "purity"
  colnames(contams) <- paste(colnames(contams), "c", sep = "_")
  
  data <- cbind(data, contams)
  
  relationshipFormula <- paste0(
    marker, "~", to 
  )
  
  modelData <- data[data[,cellType] == from, ]
  model <- lm(formula(relationshipFormula), modelData)
  
  #print(summary(model))
  
  data[data[,cellType] == from, "fittedValues"] <- predict(
    model,
    data[data$cellType == from, ]
  )
  
  
  if (plotModelFit == TRUE) {
    data[data$cellType == from, marker] <- predict(
      model, data[data$cellType == from, ]
    )
  }
  
  var <- "CD8.T_c"
  status <- "status"
  
  g1 <- ggplot2::ggplot() +
     ggplot2::stat_density_2d(
      data = data[data$cellType == to, ],
      ggplot2::aes(
        x = x, y = y, fill = ..density..
      ),
      geom = "raster",
      contour = FALSE
    )    +
    ggplot2::geom_point(data = data[data$cellType == to, ], aes(x,y), size = size/2, colour = "darkblue")+
    ggplot2::scale_fill_distiller(palette = "Blues", direction = 1) +
    ggplot2::geom_point(
      data = data[data$cellType == from, ],
      ggplot2::aes_string(
        x = "x", y = "y",
        colour = marker
      ), 
      size = size,
      shape = shape
    ) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(paste(
      "Cell Points:", from, ",",
      "Cell Density:", to, ",",
      "Model Fit:", plotModelFit
    )) +
    ggplot2::facet_wrap(~imageID, scales = "free") +
    # ggplot2::scale_color_viridis_c(option = "plasma") +
    scale_colour_gradientn(colours = rep(c("black","darkred", "red", "orange","yellow"),c(1,3,3,3,3)))
  
  
  g4 <- ggplot2::ggplot() +
    ggplot2::stat_density_2d(
      data = data[data$cellType == to, ],
      ggplot2::aes(
        x = x, y = y, fill = ..density..
      ),
      geom = "raster",
      contour = FALSE
    ) +
    ggplot2::geom_point(data = data[data$cellType == to, ], aes(x,y), size = size/2, colour = "darkblue")+
    ggplot2::scale_fill_distiller(palette = "Blues", direction = 1) +
    ggplot2::geom_point(
      data = data[data$cellType == from, ],
      ggplot2::aes_string(
        x = "x", y = "y",
        colour = var
      ), 
      size = size,
      shape = shape
    ) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(paste(
      "Cell Points:", from, ",",
      "Cell Density:", to, ",",
      "Model Fit:", plotModelFit
    )) +
    ggplot2::facet_wrap(~imageID, scales = "free") +
    # scale_colour_gradientn(colours = rep(c("black","darkred", "red", "orange","yellow"),c(1,3,3,3,3)))
    ggplot2::scale_color_viridis_c(option = "rocket")
  
  
  g2 <- data %>%
    dplyr::filter(cellType == from) |>
    ggplot2::ggplot(
      ggplot2::aes_string(
        x = to,
        y = marker,
        color = var
      )
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_viridis_c(option = "rocket") +
    ggplot2::geom_smooth(method = lm, formula = y ~ x) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("State Change Scatter Plot") +
    ggplot2::ylab(paste(marker, "expression")) +
    ggplot2::xlab(paste(from, " ", type, " to ", to)) +
    ggplot2::ylim(-1, NA)
  
  
  g3 <- data %>%
    dplyr::filter(cellType == from) |>
    ggplot2::ggplot(
      ggplot2::aes_string(
        x = to,
        y = var,
        color = marker
      )
    ) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = lm, formula = y ~ x) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("State Change Scatter Plot") +
    ggplot2::ylab(paste(marker, "expression")) +
    ggplot2::xlab(paste(from, " ", type, " to ", to)) +
    ggplot2::ylim(0, 1)
  
  
  g5 <- data %>%
    dplyr::filter(cellType == from) |>
    ggplot2::ggplot(
      ggplot2::aes_string(
        x = var,
        y = marker
      )
    ) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = lm, formula = y ~ x) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("State Change Scatter Plot") +
    ggplot2::ylab(paste(marker, "expression")) +
    ggplot2::xlab("purity") +
    ggplot2::ylim(-1, NA)

  g7 <- data %>%
    dplyr::filter(cellType == from) |>
      mutate(status = case_when(
      CD8A > 0.5 ~ "Positive",
      TRUE ~ "Negative"
    )) |>
    ggplot2::ggplot(
      ggplot2::aes_string(
        x = status,
        y = var
      )
    ) +
    ggplot2::geom_boxplot() +
    # ggplot2::geom_smooth(method = lm, formula = y ~ x) +
    ggplot2::theme_classic()
    # ggplot2::ggtitle("State Change Scatter Plot") +
    # ggplot2::ylab(paste(marker, "expression")) +
    # ggplot2::xlab("purity") +
    # ggplot2::ylim(-1, NA)


  
  g6 <- grid.arrange(g1 + xlim(0, 1500) + ylim(4000, NA), 
                     g4 + xlim(0, 1500) + ylim(4000, NA), 
                     ncol = 2)
  
  # g3 <- data |>
  #   dplyr::filter(cellType == from) |>
  #   ggplot2::ggplot(ggplot2::aes_string(
  #     x = marker,
  #     y = "fittedValues"
  #   )) +
  #   ggplot2::geom_point() +
  #   ggplot2::geom_smooth(method = lm) +
  #   ggplot2::theme_classic() +
  #   ggplot2::xlab("True Values") +
  #   ggplot2::ylab("Fitted Values") +
  #   ggplot2::ggtitle("Predicted vs Real Values")
  
  # g4 <- ggplot2::autoplot(model) + ggplot2::theme_classic()
  
  if (interactive == TRUE) {
    g1 <- plotly::ggplotly(g1)
    g2 <- plotly::ggplotly(g2)
    g3 <- plotly::ggplotly(g3)
  }
  list(image = g1, scatter = g2, contam_image = g4, marker = g5, compare = g6, boxplot = g7)
  # list(g1, g2, g3, g4)
}




#' Extract the average expression for all markers for each cell type in each 
#' region defined by lisaClust
#'
#' Takes a SingleCellExperiment and outputs a dataframe in a convenient format for 
#' cross validation
#' 
#'
#' @param data
#'   A SingleCellExperiment object with intensities data in the assays slot and
#'   regions information in colData generated by lisaClust.
#' @param imageID
#'   The colData column that stores the image IDs.
#' @param cellType
#'   The colData column that store the cell types.
#' @param region 
#'   The colData column that stores the regions.
#' @param markers 
#'   A string vector of markers that proxy a cell's state. If NULL, all markers 
#'   will be used.
#' @param assay Which assay do you want to use for the expression data.
#' @param replaceVal A value to replace missing values with.
#'
#' @examples
#' data(kerenSCE)
#' 
#' kerenSCE <- kerenSCE[,kerenSCE$imageID %in% c("5","6")]
#' 
#' regionSCE <- lisaClust::lisaClust(kerenSCE, k = 5)
#' 
#' lisaClustOutput <- getMarkerMeans(regionSCE)
#' @export
#' @rdname getMarkerMeans
#' @importFrom dplyr left_join group_by
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom SummarizedExperiment colData assay
#' @importFrom tibble column_to_rownames
getMarkerMeans <- function(data,
                       imageID = NULL,
                       cellType = NULL,
                       region = NULL,
                       markers = NULL,
                       assay = 1,
                       replaceVal = 0) {
  
  
  
  if(is.null(markers)) {
    markers <- rownames(data)
  }
  
  if(!is.null(imageID)) {
    if(!imageID%in%colnames(colData(data)))stop("Your imageID is not in colData")
    data$imageID <- colData(data)[,imageID]
    imageID <- "imageID"
  }
  
  if(!is.null(cellType)) {
    if(!cellType%in%colnames(colData(data)))stop("Your cellType is not in colData")
    data$cellType <- colData(data)[,cellType]
    cellType <- "cellType"
  }
  
  if(!is.null(region)) {
    if(!region%in%colnames(colData(data)))stop("Your region is not in colData")
    data$region <- colData(data)[,region]
    region <- "region"
  }
  
  
  
  use <- c(imageID, cellType, region)
  
  df <- data.frame(colData(data)[,use, drop = FALSE], t(assay(data, assay)))
  df <- tidyr::pivot_longer(df, -any_of(use), names_to = "markers") 
  df <- dplyr::group_by(df, across(-value))
  
  if(!is.null(imageID)){
  m <-  tidyr::pivot_wider(df, names_from = c(markers, cellType, region), values_from = value, values_fn = mean, names_sep = "__", values_fill =  replaceVal)
  m <- column_to_rownames(m, "imageID")
  m <- as.data.frame(m)
  }
  
  if(is.null(imageID)){
    m <-  tidyr::pivot_wider(df, names_from = markers, values_from = value, values_fn = mean, names_sep = "__", values_fill =  replaceVal)
    m <- as.data.frame(m)
  }
  m
}


