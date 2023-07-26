#' @noRd
#'
#' @import tidyverse
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
preProcessing <- function(SCE) {
  
  
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
  
  
  intensitiesData <- data.frame(t(assay(SCE, "intensities")))
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
    colData(SCE)[, x] <- colData(SCE)[, x] %>% as.character
  }

  return(SCE)
}


#' Calculate pairwise distance between cell types
#'
#' Calculates the euclidean distance from each cell to the nearest cell of each
#' type for a single image
#'
#' @param singleCellData the single cell data of interest
#' @param maxRS Maximum distance between pairs of points to be counted as close
#'   pairs.
#'
#' @export
#' @rdname distanceCalculator
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr
#'   select mutate inner_join arrange group_by slice rename relocate full_join
#' @importFrom spatstat.geom owin ppp closepairs
#' @importFrom tidyr pivot_wider
distanceCalculator <- function(singleCellData, maxRS = 200) {
  singleCellData <- singleCellData %>%
    tibble::rownames_to_column("i") %>%
    dplyr::select(-i) %>%
    tibble::rownames_to_column("cellIndex") %>%
    dplyr::mutate(cellIndex = as.numeric(cellIndex))
  
  
  ow <- spatstat.geom::owin(
    xrange = range(singleCellData$x),
    yrange = range(singleCellData$y)
  )
  pppData <- spatstat.geom::ppp(
    x = singleCellData$x,
    y = singleCellData$y,
    window = ow,
    marks = singleCellData$cellType
  )
  
  closePairData <- spatstat.geom::closepairs(pppData, rmax = maxRS)
  distanceData <- data.frame(
    cellIndexA = closePairData$i,
    cellIndexB = closePairData$j,
    d = closePairData$d
  )
  
  
  cellAInformation <- singleCellData %>%
    dplyr::select(
      cellIndexA = cellIndex,
      cellTypeA = cellType,
      cellIDA = cellID
    )
  cellBInformation <- singleCellData %>%
    dplyr::select(cellIndexB = cellIndex, cellTypeB = cellType)
  
  processedDistanceData <- distanceData %>%
    dplyr::inner_join(cellAInformation, by = "cellIndexA") %>%
    dplyr::inner_join(cellBInformation, by = "cellIndexB") %>%
    dplyr::arrange(d) %>%
    dplyr::group_by(cellIndexA, cellTypeA, cellTypeB) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::rename(cellType = "cellTypeA") %>%
    dplyr::rename(cellIndex = "cellIndexA") %>%
    dplyr::rename(cellID = "cellIDA") %>%
    dplyr::select(-cellIndexB, -cellIndex) %>%
    dplyr::mutate(cellTypeB = paste0("dist_", cellTypeB)) %>%
    tidyr::pivot_wider(names_from = cellTypeB, values_from = d) %>%
    dplyr::mutate(imageID = unique(singleCellData$imageID)) %>%
    dplyr::relocate(imageID)
  
  
  processedDistanceData <- processedDistanceData %>%
    dplyr::right_join(singleCellData, by = c("imageID", "cellID", "cellType"))
}


#' Wrapper to calculate pairwise distance between cell types by image
#'
#' Calculates the euclidean distance from each cell to the nearest cell of each
#' type
#'
#' @param singleCellData
#'   A dataframe with a cellType column as well as x and y spatial coordinates.
#'   The dataframe must contain a imageID column and cellID (unique cell
#'   identifier's) column as well
#' @param Rs
#'   Radius to calculate pairwise distances between cells (can be a numeric or
#'   vector of radii)
#' @param whichCellTypes
#'   Character vector specifying what cell types to include in the calculation.
#'   If the argument is non-null, then at least two cell types need to be
#'   specified
#' @param nCores Number of cores for parallel processing
#'
#' @examples
#' library(dplyr)
#' data("headSCE")
#'
#' singleCellDataDistances <- getDistances(headSCE,
#'   nCores = 1,
#'   Rs = c(200),
#'   whichCellTypes = c("MC2", "SC7")
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
getDistances <- function(singleCellData,
                         Rs = c(200),
                         whichCellTypes = NULL,
                         nCores = 1) {
  x <- runif(1)
  BPPARAM <- .generateBPParam(cores = nCores)
  
  markersToUse <- rownames(singleCellData)
  
  metadata_name <- paste0("Rs", Rs, "")
  
  if(!metadata_name %in% names(metadata(singleCellData))) {
    
    if(!all(markersToUse %in% colnames(colData(singleCellData)))) {
      
      singleCellDataClean <- singleCellData %>%
        preProcessing()
      
      singleCellData <- as.data.frame(colData(singleCellDataClean))
      
    }
    
  } else {
    singleCellDataClean <- singleCellData
    singleCellData <- metadata(singleCellDataClean)[[metadata_name]]
  }
  
  if (!is.null(whichCellTypes)) {
    if (length(whichCellTypes) >= 2) {
      singleCellData <- singleCellData %>%
        dplyr::filter(cellType %in% whichCellTypes)
    } else {
      warnings(
        "The argument whichCellTypes needs to contain at least cellType",
        " names in the vector"
      )
    }
  }
  
  singleCellDataDistances <- singleCellData %>%
    split(~imageID) %>%
    BiocParallel::bplapply(distanceCalculator,
                           maxRS = max(Rs),
                           BPPARAM = BPPARAM
    ) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(dplyr::across(where(is.numeric), function(x) ifelse(is.infinite(x), NA, x))) %>%
    lapply(Rs, function(x, rmax) {
      x %>%
        dplyr::mutate(
          dplyr::across(dplyr::contains("dist_"),
                        function(x) ifelse(x <= rmax, x, NA)
          )
        ) %>%
        dplyr::rename_with(
          function(x) {
            stringr::str_replace(x, "dist_", paste0("dist", rmax, "_"))
          },
          dplyr::starts_with("dist_")
        )
    },
    x = .
    ) %>%
    purrr::reduce(full_join)
  
  metadata_name <- paste0("Rs", Rs, "")
  metadata(singleCellDataClean)[[metadata_name]] <- singleCellDataDistances
  # # Identify overlapping column names
  # metadata(singleCellDataClean) <- append(metadata(singleCellDataClean), list(singleCellDataDistances))
  # overlap_cols <- intersect(names(singleCellDataDistances), colnames(colData(singleCellDataClean)))
  # 
  # # Remove overlapping columns from the new data frame
  # singleCellDataDistances <- singleCellDataDistances[, !names(singleCellDataDistances) %in% overlap_cols]
  # 
  # # Append singleCellDataDistances to colData slot in SCE
  # colData(singleCellDataClean) <- cbind(colData(singleCellDataClean), singleCellDataDistances)
  # 
  return(singleCellDataClean)
}



#' Wrapper to calculate imhomogenous K function between a cell and surrounding
#' types on each image
#'
#' Calculate the imhomogenous K function (a measure of cell type abundance) for
#' each cell to other cell types
#'
#' @param singleCellData
#'   A dataframe with a cellType column as well as x and y spatial coordinates.
#'   The dataframe must contain a imageID column and cellID (unique cell
#'   identifier's) column as well
#' @param Rs
#'   Radius to include in that calculation of pairwise abundance (K-function)
#'   between cells (can be a numeric or vector of radii)
#' @param whichCellTypes
#'   Character vector specifying what cell types to include in the calculation.
#'   If the argument is non-null, then at least two celltypes must be specified.
#' @param nCores Number of cores for parallel processing
#'
#' @examples
#' library(dplyr)
#' data("headSCE")
#' 
#' singleCellDataCounts <- getAbundances(headSCE,
#'   nCores = 1,
#'   Rs = c(200),
#'   whichCellTypes = c("MC2", "SC7")
#' )
#'
#' @export
#' @rdname getAbundances
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr across bind_rows inner_join relocate contains mutate
#' @importFrom stringr word
#' @importFrom SummarizedExperiment colData assayNames
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom magrittr %>%
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
getAbundances <- function(singleCellData,
                          Rs = c(200),
                          whichCellTypes = NULL,
                          nCores = 1) {
  
  x <- runif(1)
  BPPARAM <- .generateBPParam(cores = nCores)
  
  markersToUse <- rownames(singleCellData)
  
  metadata_name <- paste0("Rs", Rs, "")
  
  if(!metadata_name %in% names(metadata(singleCellData))) {
    
    if(!all(markersToUse %in% colnames(colData(singleCellData)))) {
      
      singleCellDataClean <- singleCellData %>%
        preProcessing()
      
      singleCellData <- as.data.frame(colData(singleCellDataClean))
    
    }
    
  } else {
    singleCellDataClean <- singleCellData
    singleCellData <- metadata(singleCellDataClean)[[metadata_name]]
  }
  
  if (!is.null(whichCellTypes)) {
    if (length(whichCellTypes) >= 2) {
      singleCellData <- singleCellData %>%
        dplyr::filter(cellType %in% whichCellTypes)
    } else {
      warnings(
        "The argument whichCellTypes needs to contain at least cellType",
        " names in the vector"
      )
    }
  }
  
  # Sourish: Make sure to have check to ensure rownames are not purely numbers
  
  singleCellDataK <- singleCellData %>%
    split(~imageID) %>%
    BiocParallel::bplapply(
      lisaClust:::inhomLocalK,
      Rs = Rs,
      BPPARAM = BPPARAM
    ) %>%
    lapply(as.data.frame)
  
  
  # Correcting column names when Rs greater than max Rs of image
  # Nick: this is a nested nightmare 😨
  # Sourish: Haha, I think this is a bug from the output of inhomLocalK from the lisaClust package that I tried to fix
  correctedRadius <- singleCellDataK %>%
    lapply(
      function(x) {
        unlist(
          lapply(
            colnames(x),
            function(x) {
              min(Rs[Rs >= as.numeric(stringr::word(x, 1, sep = "_"))])
            }
          )
        )
      }
    )
  
  cellTypeName <- singleCellDataK %>%
    lapply(function(x) word(colnames(x), 2, -1, sep = "_"))
  correctedColumnNames <- correctedRadius %>%
    mapply(
      function(newRadius, cellTypeName) {
        mapply(function(a, b) paste0(a, "_", b), newRadius, cellTypeName)
      },
      newRadius = .,
      cellTypeName = cellTypeName,
      SIMPLIFY = FALSE
    )
  
  singleCellDataK <- singleCellDataK %>%
    mapply(
      function(x, newNames) x %>% purrr::set_names(newNames),
      x = ., newNames = correctedColumnNames, SIMPLIFY = FALSE
    ) 
  singleCellDataK <- singleCellDataK %>%
    mapply(
      function(x, imageID) {
        x %>%
          dplyr::rename_all(~ paste0("abundance", .x)) %>%
          dplyr::mutate(imageID = imageID)
      },
      x = ., imageID = names(.), SIMPLIFY = FALSE
    ) 
  singleCellDataK <- singleCellDataK %>%
    dplyr::bind_rows() %>%
    tibble::rownames_to_column("cellID") %>%
    dplyr::mutate(dplyr::across(where(is.numeric), function(x) ifelse(is.na(x), 0, x)))
  
  singleCellDataK <- singleCellDataK %>%
    dplyr::right_join(singleCellData, by = c("imageID", "cellID")) %>%
    dplyr::relocate(imageID)
  
  metadata(singleCellDataClean)[[metadata_name]] <- singleCellDataK
  return(singleCellDataClean)
}


#' Calculate the level of marker contamination of each cell
#'
#' Calculates contamination scores using a random forest classification
#'
#' @param singleCellData
#'   A dataframe with a cellType column as well as marker intensity information
#'   corresponding to each cell. The dataframe must contain a imageID column.
#' @param Rs
#'   Radius to calculate pairwise distances between cells (can be a numeric or
#'   vector of radii)
#' @param seed
#'   A numeric value to allow for replicable results from the random forest
#'   classification algorithm
#' @param num.trees Number of trees to be used in the random forest classifier
#' @param verbose
#'   A logical indicating whether information about the final random forest
#'   model should be outputted.
#' @param missingReplacement
#'   A default value to replace missing marker intensities for classification.
#'
#' @examples
#' library(dplyr)
#' data("headSCE")
#'
#' singleCellDataDistancesContam <- calcContamination(
#'   headSCE,
#'   Rs = c(200)
#' )
#'
#' @export
#' @rdname calcContamination
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate select bind_rows inner_join across
#' @importFrom SummarizedExperiment colData assayNames
#' @importFrom ranger ranger
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom stringr str_replace
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
calcContamination <- function(singleCellData,
                              Rs = c(200),
                              seed = 2022,
                              num.trees = 100,
                              verbose = FALSE,
                              missingReplacement = 0) {
  # if ("SingleCellExperiment" %in% class(singleCellData)) {
  #   singleCellDataNew <- data.frame(
  #     SummarizedExperiment::colData(singleCellData)
  #   ) %>%
  #     dplyr::mutate(
  #       imageID = as.character(imageID), cellType = as.character(cellType)
  #     )
  # 
  #   if ("intensities" %in% assayNames(singleCellData)) {
  #     singleCellDataNew <- singleCellDataNew %>%
  #       cbind(t(SummarizedExperiment::assay(singleCellData, "intensities")))
  #     if (is.null(markers)) {
  #       markers <- rownames(
  #         SummarizedExperiment::assay(singleCellData, "intensities")
  #       )
  #     }
  #   }
  # 
  # 
  #   singleCellData <- singleCellDataNew
  # }
  
  
  markers <- rownames(singleCellData)
  
  metadata_name <- paste0("Rs", Rs, "")
  
  if(!all(markers %in% colnames(colData(singleCellData)))) {
    
    singleCellDataClean <- singleCellData %>%
      preProcessing()
    
    singleCellData <- as.data.frame(colData(singleCellDataClean))
    
  } else {
    
    singleCellDataClean <- singleCellData
    singleCellData <- as.data.frame(colData(singleCellDataClean))  
  
  }
  
  
  
  rfData <- singleCellData %>%
    dplyr::select(cellType, all_of(markers)) %>%
    dplyr::mutate(dplyr::across(markers, function(x) ifelse(is.nan(x) | is.na(x), 0, x)))
  
  set.seed(seed)
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


  if(!metadata_name %in% names(metadata(singleCellDataClean))) {
    
    distData <- colData(singleCellDataClean) 
    
    distData <- distData %>% as.data.frame() %>%
      dplyr::left_join(rfData)
    
    metadata(singleCellDataClean)[[metadata_name]] <- distData
    
  } else {
    
    distData <- metadata(singleCellDataClean)[[metadata_name]]
    
    distData <- distData %>% as.data.frame() %>%
      dplyr::left_join(rfData)
    
    metadata(singleCellDataClean)[[metadata_name]] <- distData
    
  }
  
  return(singleCellDataClean)
}







### Section 1: LM Calculations ##############################################################


#' First layer wrapper function to build linear models measuring state changes
#'
#' Builds linear models measuring marker based state changes in a cell type
#' based of the proximity or abundance of another cell type. The function
#' provides the option to build robust and mixed linear model variants
#'
#' @param singleCellData
#'   A dataframe with a imageID, cellType, and marker intensity column along
#'   with covariates (e.g. distance or abundance of the nearest cell type) to
#'   model cell state changes
#' @param Rs
#'   Radius to calculate pairwise distances between cells (can be a numeric or
#'   vector of radii)
#' @param typeAll
#'   A prefix that appears on the column names of all cell state modelling
#'   covariates. The default value is "dist"
#' @param covariates
#'   A list of additional covariates to be included in the model being built.
#' @param condition
#'   A list of additional conditions to be included in the model being built.
#' @param method
#'  The type of linear model to be built. Current options include "lm", "rlm".
#' @param isMixed
#'   A logical indicating if a mixed effects model should be built - this will
#'   build models for relationships on a global basis rather than on a image by
#'   image level
#' @param randomIntercepts
#'   A string list of column names in the dataframe that should correspond to
#'   the random intercepts
#' @param cellTypesToModel
#'   A string vector of the names of cell types to model cell state changes in.
#'   The default argument is NULL which models are cell types
#' @param verbose A logical indicating if messages should be printed
#' @param timeout
#'   A maximum time allowed to build each model. Setting this may be important
#'   when building rlm mixed linear models
#' @param nCores Number of cores for parallel processing
#'
#' @examples
#' library(dplyr)
#' data("headSCE")
#'
#' singleCellDataDistances <- getDistances(headSCE,
#'   nCores = 1,
#'   Rs = c(200),
#'   whichCellTypes = c("MC2", "SC7")
#' )
#'
#' imageModels <- getStateChanges(
#'   singleCellData = singleCellDataDistances,
#'   Rs = c(200),
#'   typeAll = c("dist200"),
#'   cellTypesToModel = "MC2",
#'   nCores = 1
#' )
#'
#' @export
#' @rdname getStateChanges
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom dplyr
#'   arrange group_by  summarise_at mutate bind_rows left_join filter
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
getStateChanges <- function(singleCellData,
                            Rs,
                            typeAll = c("dist"),
                            covariates = NULL,
                            condition = NULL,
                            method = "lm",
                            isMixed = FALSE,
                            randomIntercepts = c("imageID"),
                            cellTypesToModel = NULL,
                            verbose = FALSE,
                            timeout = 10,
                            nCores = 1) {
  
  markers = rownames(singleCellData)
  
  metadata_name <- paste0("Rs", Rs, "")
  
  SCE <- singleCellData 
  singleCellData <- metadata(SCE)[[metadata_name]]
  
  typeVector <- singleCellData %>%
    dplyr::select(contains(typeAll)) %>%
    colnames(.) %>%
    stringr::str_split("_") %>%
    lapply(function(x) x[[1]]) %>%
    unlist() %>%
    unique() %>%
    paste0("_")
  
  if (!is.null(cellTypesToModel)) {
    singleCellData <- singleCellData %>%
      dplyr::filter(cellType %in% cellTypesToModel)
  }
  
  
  allModels <- typeVector %>%
    mapply(calculateStateModels, type = ., MoreArgs = list(
      singleCellData = singleCellData,
      markers = markers,
      covariates = covariates,
      condition = condition,
      method = method,
      isMixed = isMixed,
      randomIntercepts = randomIntercepts,
      verbose = verbose,
      timeout = timeout,
      nCores = nCores
    ), SIMPLIFY = FALSE) %>%
    dplyr::bind_rows()
}




#' Second layer wrapper function to build linear models measuring state changes
#'
#' Calculates linear models For pairwise interactions with options of mixed
#' models and robust regression with a single "type" of modelling covariate
#'
#' @noRd
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom dplyr
#'   arrange group_by  summarise_at across mutate bind_rows left_join
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
calculateStateModels <- function(singleCellData,
                                 markers,
                                 type = "dist200_",
                                 covariates = NULL,
                                 condition = NULL,
                                 method = "lm",
                                 isMixed = FALSE,
                                 randomIntercepts = c("imageID"),
                                 verbose = TRUE,
                                 timeout = 10,
                                 nCores = 1) {
  
  x <- runif(1)
  BPPARAM <- .generateBPParam(cores = nCores)
  
  if (isMixed == TRUE) {
    splitData <- split(singleCellData, ~cellType)
  } else {
    splitData <- split(singleCellData, ~imageID)
  }
  
  CellInteractionModels <- splitData %>%
    BiocParallel::bplapply(buildModelsByCellType,
                           markers = markers,
                           covariates = covariates,
                           condition = condition,
                           type = type,
                           method = method,
                           isMixed = isMixed,
                           randomIntercepts = randomIntercepts,
                           verbose = verbose,
                           timeout = timeout,
                           BPPARAM = BPPARAM
    ) %>%
    bind_rows()
  
  
  CellInteractionModels <- CellInteractionModels %>%
    dplyr::arrange(pValue) %>%
    mutate(type = type) %>%
    mutate(covariates = paste(covariates, collapse = " + ")) %>%
    mutate(covariates = ifelse(is.null(covariates), "None", covariates))
  
  relativeExpressionData <- singleCellData %>%
    dplyr::group_by(independent = cellType) %>%
    dplyr::summarise_at(markers, mean, na.rm = TRUE) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), heatmaply::normalize)) %>%
    tidyr::gather(
      -independent,
      key = "dependent", value = "relativeExpression"
    ) %>%
    dplyr::mutate(interactingCell = independent) %>%
    dplyr::mutate(independent = paste0(type, independent))
  
  CellInteractionModels <- CellInteractionModels %>%
    dplyr::mutate(interactingCell = str_replace(independent, type, "")) %>%
    dplyr::left_join(relativeExpressionData) %>%
    relocate(cellType, independent, dependent)
  
  if (isMixed == FALSE) {
    CellInteractionModels <- CellInteractionModels %>%
      relocate(imageID)
  }
  
  CellInteractionModels
}


#' Third layer wrapper function to initiate building linear models measuring
#' state changes on a cell type by cell type basis
#'
#' @noRd
#'
#' @importFrom dplyr group_by  group_modify ungroup
#' @importFrom magrittr %>%
buildModelsByCellType <- function(subsettedSingleCellData,
                                  markers,
                                  type,
                                  covariates,
                                  condition,
                                  method,
                                  isMixed,
                                  randomIntercepts,
                                  verbose,
                                  timeout) {
  finalModelsOutputs <- subsettedSingleCellData %>%
    dplyr::group_by(cellType) %>%
    dplyr::group_modify(
      ~ .x %>% modelsPerCellType(
        markers, type, covariates, condition, method, isMixed,
        randomIntercepts, verbose, timeout
      )
    ) %>%
    dplyr::ungroup()
}




#' Fourth layer wrapper function to initiate creating model formulas and
#' building the state change detection linear models
#'
#' @noRd
#'
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_detect str_replace
#' @importFrom magrittr %>%
modelsPerCellType <- function(cellTypeSplitData,
                              markers,
                              type,
                              covariates,
                              condition,
                              method,
                              isMixed,
                              randomIntercepts,
                              verbose,
                              timeout) {
  if (verbose == TRUE) {
    print(unique(cellTypeSplitData$imageID))
    print(unique(cellTypeSplitData$cellType))
  }
  
  cells <- colnames(cellTypeSplitData)[
    Reduce(
      "|", mapply(
        function(x, type) stringr::str_detect(colnames(x), type),
        type,
        MoreArgs = list(x = cellTypeSplitData),
        SIMPLIFY = FALSE
      )
    )
  ]
  
  # Univariate Models
  if (length(type) == 1) {
    formulas <- lapply(markers, function(x, c) paste(x, "~", c), c = cells) %>%
      unlist()
    if (!is.null(condition)) {
      formulas <- formulas %>% lapply(function(x, c) paste(x, "*", c), c = condition) %>% unlist()  
    }
    
    # Interaction Models - two types
  } else {
    for (i in type) {
      cells <- stringr::str_replace(cells, i, "")
    }
    cells <- unique(cells)
    cells <- lapply(cells, function(x) unlist(lapply(type, paste0, x))) %>%
      lapply(
        function(x) {
          paste0(
            c(
              paste0(x, collapse = " + "),
              paste0(x, collapse = ":")
            ),
            collapse = " + "
          )
        }
      ) %>%
      unlist()
    formulas <- lapply(markers, function(x, c) paste(x, "~", c), c = cells) %>%
      unlist()
  }
  
  modelOutputs <- formulas %>%
    mapply(fitStateModels,
           f = .,
           MoreArgs = list(
             x = cellTypeSplitData,
             covariates = covariates,
             method = method,
             isMixed = isMixed,
             randomIntercepts = randomIntercepts,
             timeout = timeout
           ),
           SIMPLIFY = FALSE
    ) %>%
    dplyr::bind_rows()
  
  modelOutputs
}





#' Fifth layer wrapper that builds the cell state models
#'
#' @noRd
#'
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_detect str_replace str_split
#' @importFrom R.utils withTimeout
#' @importFrom robustlmm rlmer
#' @importFrom parameters p_value
#' @importFrom lmerTest lmer
#' @importFrom MuMIn r.squaredGLMM
#' @importFrom MASS rlm
#' @importFrom sfsmisc f.robftest
#' @importFrom performance check_singularity
#' @importFrom magrittr %>%
fitStateModels <- function(x,
                           f,
                           covariates,
                           method,
                           isMixed,
                           randomIntercepts,
                           timeout) {
  dependent <- stringr::str_split(f, " ~ ")[[1]][1]
  independent <- stringr::str_split(f, " ~ ")[[1]][2]
  
  independentCheck <- stringr::str_split(independent, "\\:") %>%
    unlist() %>%
    stringr::str_split(" \\+ ") %>%
    unlist() %>%
    unique()
  
  f <- paste(c(f, covariates), collapse = " + ")
  independentSplit <- unlist(stringr::str_split(independent, " \\* "))
  
  outputs <- try(
    {
      if (isMixed == TRUE) {
        randomInterceptTerms <- randomIntercepts %>%
          lapply(function(x) paste0("(1|", x, ")")) %>%
          unlist() %>%
          paste0(collapse = " + ")
        
        if (method == "rlm") {
          model <- R.utils::withTimeout(
            {
              robustlmm::rlmer(formula(paste(f, "+", randomInterceptTerms)),
                               method = "DASvar",
                               data = x
              )
            },
            timeout = timeout,
            onTimeout = "silent"
          )
          coefs <- coef(summary(model))
          coefs <- coefs[, colnames(coefs) != "df"]
          coefs <- cbind(coefs, parameters::p_value(model)$p)
          modelSummary <- summary(model)
          modelSummary$r.squared <- NA
        } else {
          model <- R.utils::withTimeout(
            {
              lmerTest::lmer(formula(paste(f, "+ (1|imageID)")), data = x)
            },
            timeout = timeout,
            onTimeout = "silent"
          )
          coefs <- coef(summary(model))
          coefs <- coefs[, colnames(coefs) != "df"]
          modelSummary <- summary(model)
          modelSummary$r.squared <- MuMIn::r.squaredGLMM(model)[1]
        }
      } else {
        if (method == "rlm") {
          model <- R.utils::withTimeout(
            {
              MASS::rlm(formula(f), data = x)
            },
            timeout = timeout,
            onTimeout = "silent"
          )
          coefs <- coef(summary(model))
          pValues <- unlist(
            lapply(
              independentSplit,
              function(x) sfsmisc::f.robftest(model, var = x)$p.value
            )
          )
          names(pValues) <- independentSplit
          coefs <- cbind(coefs, data.frame(pValue = NA))
          coefs[names(pValues), "pValue"] <- pValues
          modelSummary <- summary(model)
        } else {
          model <- R.utils::withTimeout(
            lm(formula(f), data = x),
            timeout = timeout, onTimeout = "silent"
          )
          coefs <- coef(summary(model))
          modelSummary <- summary(model)
        }
      }
      
      beta <- coefs[independentSplit[1], 1]
      tValue <- coefs[independentSplit[1], 3]
      pValue <- coefs[independentSplit[1], 4]
      rValue <- modelSummary$r.squared
      sampleSize <- length(modelSummary$residuals)
      outputs <- data.frame(
        beta = beta,
        tValue = tValue,
        pValue = pValue,
        independent = paste0(independentSplit, collapse = ", "),
        dependent = dependent,
        rValue = rValue,
        sampleSize = sampleSize,
        formula = f,
        isSingular = performance::check_singularity(model)
      )
      if (isMixed == FALSE) {
        outputs <- outputs %>%
          mutate(imageID = unique(x$imageID))
      }
      outputs
    },
    silent = TRUE
  )
  
  if (any(class(outputs) == "try-error")) {
    outputs <- data.frame(
      independent = paste0(independentSplit, collapse = ""),
      dependent = dependent,
      formula = f
    )
    
    if (isMixed == FALSE) {
      outputs <- outputs %>%
        mutate(imageID = unique(x$imageID))
    }
  }
  
  outputs
}

### Section 2:Fast Version LM ##############################################################


#' Wrapper function to quickly build ordinary linear models measuring state
#' changes on a image basis
#'
#' Builds linear models measuring marker based state changes in a cell type
#' based of the proximity or abundance of another cell type. The function only
#' provides the option to build OLS models
#'
#' @param singleCellData
#'   A dataframe with a imageID, cellType, and marker intensity column along
#'   with covariates (e.g. distance or abundance of the nearest cell type) to
#'   model cell state changes
#' @param Rs Radius to calculate pairwise distances between cells (can be a numeric or
#'   vector of radii)
#' @param type
#'   A prefix that appears on the column names of all cell state modelling
#'   covariates. The default value is "dist"
#' @param covariates
#'   A list of additional covariates to be included in the model being built
#' @param removeColsThresh
#'   A threshold in which a modelling covariate column will be excluded in a
#'   image model if missingness exceeds a certain proportion
#' @param cellTypesToModel
#'   A string vector of the names of cell types to model cell state changes in.
#'   The default argument is NULL which models are cell types
#' @param nCores Number of cores for parallel processing
#'
#'
#' @examples
#' library(dplyr)
#' data("headSCE")
#'
#' singleCellDataDistances <- getDistances(headSCE,
#'   nCores = 1,
#'   Rs = c(200),
#'   whichCellTypes = c("MC2", "SC7")
#' )
#'
#' imageModelsFast <- getStateChangesFast(
#'   singleCellData = singleCellDataDistances,
#'   Rs = c(200),
#'   type = c("dist200"),
#'   nCores = 1
#' )
#' @export
#' @rdname getStateChangesFast
#' @importFrom dplyr bind_rows filter
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom magrittr %>%
getStateChangesFast <- function(singleCellData,
                                Rs,
                                type = "dist",
                                covariates = NULL,
                                removeColsThresh = 0.1,
                                cellTypesToModel = NULL,
                                nCores = 1) {
  
  x <- runif(1)
  BPPARAM <- .generateBPParam(cores = nCores)
  
  markers = rownames(singleCellData)
  
  metadata_name <- paste0("Rs", Rs, "")
  
  SCE <- singleCellData 
  singleCellData <- metadata(SCE)[[metadata_name]]
  
  if (!is.null(cellTypesToModel)) {
    singleCellData <- singleCellData %>%
      dplyr::filter(cellType %in% cellTypesToModel)
  }
  
  imageModels <- BiocParallel::bplapply(
    split(singleCellData, ~ imageID + cellType),
    calculateStateModelsFast,
    type = type,
    markers = markers,
    covariates = covariates,
    removeColsThresh = removeColsThresh,
    BPPARAM = BPPARAM
  )
  
  imageModels <- imageModels[unlist(lapply(
    imageModels,
    function(x) class(x) != "try-error"
  ))] %>%
    dplyr::bind_rows()
  
  imageModels
}




#' Fast function build state change detection OLS linear models
#'
#' @noRd
#'
#' @importFrom dplyr
#'   select starts_with all_of select_if mutate bind_rows relocate rename
#' @importFrom stringr word
#' @importFrom bigstatsr as_FBM big_univLinReg
#' @importFrom magrittr %>%
calculateStateModelsFast <- function(singleCellData,
                                     type = "dist",
                                     markers,
                                     covariates = NULL,
                                     removeColsThresh = 0.1) {
  model <- try({
    singleCellData <- singleCellData %>%
      dplyr::select(
        imageID,
        cellType,
        dplyr::starts_with(type),
        dplyr::all_of(markers),
        dplyr::all_of(covariates)
      )
    
    singleCellData <- na.omit(singleCellData[, colSums(is.na(singleCellData)) < removeColsThresh * nrow(singleCellData)])
    
    x.train.original <- singleCellData %>%
      dplyr::select(dplyr::starts_with(type))
    x.train <- bigstatsr::as_FBM(x.train.original)
    y.train <- singleCellData %>%
      dplyr::select(dplyr::all_of(markers)) %>%
      dplyr::select_if(~ length(unique(.)) > 1)
    
    y.train <- apply(y.train, MARGIN = 2, function(x) x, simplify = FALSE)
    cov.train <- singleCellData %>%
      dplyr::select(dplyr::all_of(covariates))
    
    model <- mapply(
      FUN = function(y.train) {
        bigstatsr::big_univLinReg(
          X = x.train,
          y.train = y.train,
          covar.train = as.matrix(cov.train)
        )
      },
      y.train = y.train,
      SIMPLIFY = FALSE
    ) %>%
      mapply(
        FUN = function(x, dependent) {
          x %>% dplyr::mutate(
            independent = colnames(x.train.original),
            dependent = dependent
          )
        },
        dependent = names(.),
        SIMPLIFY = FALSE
      ) %>%
      dplyr::bind_rows()
    
    model <- model %>%
      dplyr::mutate(sampleSize = nrow(x.train)) %>%
      dplyr::mutate(
        imageID = unique(singleCellData$imageID),
        cellType = unique(singleCellData$cellType),
        pValue = ifelse(
          is.na(score),
          NA,
          2 * pt(
            abs(score),
            df = sampleSize - 2 - ncol(cov.train),
            lower.tail = FALSE
          )
        )
      ) %>%
      dplyr::relocate(imageID, cellType, independent, dependent) %>%
      dplyr::mutate(
        interactingCell = word(independent, 2, -1, "_"),
        type = stringr::word(independent, 1, 1, "_")
      ) %>%
      dplyr::rename(beta = "estim", tValue = "score") %>%
      data.frame() %>%
      dplyr::mutate(covariateType = covariates)
  })
  
  model
}




### Section 3: Cross Validation ############################################################## 




#' @noRd
#' 
#' @importFrom dplyr bind_rows across mutate rename select
#' @importFrom stringr str_detect str_replace str_split
#' @importFrom tidyr pivot_wider
#' @importFrom magrittr %>%
imageModelsCVFormat <- function(imageModels,
                                values_from = "tValue",
                                removeColsThresh = 0.2,
                                missingReplacement = 0) {
  
  cvData <- imageModels %>%
    dplyr::mutate(dplyr::across(where(is.numeric), function(x) ifelse(is.finite(x), x, NA))) %>%
    dplyr::mutate(
      relationship = paste0(cellType, "_", dependent, "_", independent)
    ) %>%
    dplyr::rename(values_from = values_from) %>%
    dplyr::select(imageID, relationship, values_from) %>%
    tidyr::pivot_wider(names_from = relationship, values_from = values_from)
  
  cvData <- cvData[, colSums(is.na(cvData)) < nrow(cvData) * removeColsThresh]
  cvData[is.na(cvData)] <- missingReplacement
  
  cvData
}

#' Convert Image Model Output to Cross-Validation Format
#'
#' Takes the output from the getStateChanges function for image based state
#' models and converts it in a convenient format for cross validation
#'
#' @param imageModels
#'   A dataframe with the output from the function getStateChanges or
#'   getStateChangesFast with the argument isMixed = FALSE
#' @param values_from
#'   Column to use from imageModels for cross validation. The default is
#'   "tValue" but "beta" can be choosen as well
#' @param removeColsThresh
#'   Threshold of missingness in which a relationship will not be included as
#'   column in the cross validation ready output
#' @param missingReplacement Numeric value to replace missing values
#'
#' @examples
#' library(dplyr)
#' data("headSCE")
#'
#' singleCellDataDistances <- getDistances(headSCE,
#'   nCores = 1,
#'   Rs = c(200),
#'   whichCellTypes = c("MC2", "SC7")
#' )
#'
#' imageModels <- getStateChanges(
#'   singleCellData = singleCellDataDistances,
#'   Rs = c(200),
#'   type = c("dist200"),
#'   nCores = 1
#' )
#' crossValidationData <- listImageModelsCVFormat(imageModels,
#'   values_from = "tValue",
#'   removeColsThresh = 0.2,
#'   missingReplacement = 0
#' )
#' @export
#' @rdname imageModelsCVFormat
#' @importFrom dplyr bind_rows across mutate rename select
#' @importFrom stringr str_detect str_replace str_split
#' @importFrom tidyr pivot_wider
#' @importFrom magrittr %>%
listImageModelsCVFormat <- function(imageModels,
                                values_from = "tValue",
                                removeColsThresh = 0.2,
                                missingReplacement = 0) {
  
  if(!c("imageID") %in% colnames(imageModels)) {
    stop("The argument imageID needs to exist in the data i.e. use isMixed = FALSE")
  }                                                                             
  classificationData <- imageModels$imageID %>% unique()
  
  crossValidateInteractionsData <- imageModels %>%
    filter(sampleSize >= 10) %>%
    split(~ covariates + type) %>% bplapply(function(x) imageModelsCVFormat(x,
                                                                            values_from = values_from,
                                                                            removeColsThresh = removeColsThresh,
                                                                            missingReplacement = missingReplacement)) %>%
    mapply(function(x, dataset) x %>% mutate(dataset = dataset),
           x = .,
           dataset = names(.),
           SIMPLIFY = FALSE )
  crossValidateInteractionsData <- crossValidateInteractionsData[unlist(lapply(crossValidateInteractionsData,
                                                                              function(x) ncol(x) > 2))]
  
  modellingData <- crossValidateInteractionsData %>% 
    lapply(function(x, imageSubset) x %>%
             filter(imageID %in% imageSubset) %>% 
             column_to_rownames("imageID") %>% 
             dplyr::select(!contains(c("dataset", "type"))), 
           imageSubset = names(classificationData))
  
  return(modellingData)
}

### Section 4:Visualise Image ##############################################################

#' Visualise Cell-Cell Marker Relationships
#'
#' Helper functions to visualise OLS model fits for image based state models
#'
#' @param data
#'   A dataframe with a imageID, cellType, and marker intensity column along
#'   with covariates (e.g. distance or abundance of the nearest cell type) to
#'   model cell state changes
#' @param Rs
#'   Radius to calculate pairwise distances between cells (can be a numeric or
#'   vector of radii)
#' @param imageID
#'   Identifier name of the image in the imageID column to be visualised
#' @param mainCellType
#'   String indicating the name of the cell type (from the cellType column) whose cell state is being investigated in
#' @param interactingCellType
#' String indicating the name of the cell type (from the cellType column) who may be influencing the cell state of another cell type
#' @param depedentMarker
#'   String refering to the marker column proxying the cell state of interest
#' @param sizeVariable
#'   Aesthetic numerical variable determining the size of the displayed cells
#' @param shape
#'   Aesthetic variable determining the shape grouping of the displayed cells
#' @param interactive
#'   Logical indicating if the output visualisation should be a interactive (plotly)
#' @param plotModelFit
#'   Logical indicating if fitted values should be plotted or actual intensities for marker specified. The default is to plot actual intensities
#' @param method 
#'   The method to build the model with. Currently the only option is "lm". However, capabilities may be expanded in the future
#' @param modelType
#'   String indicating the prefix of the independent metric corresponding to the potential state change inducing cell type. For example, values could include "abundance200_" or "dist200_"
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' data("headSCE")
#'
#' singleCellDataDistances <- getDistances(headSCE,
#'   nCores = 1,
#'   Rs = c(200),
#'   whichCellTypes = c("MC2", "SC7")
#' )
#' visualiseImageRelationship(
#'   data = singleCellDataDistances,
#'   Rs = c(200),
#'   imageID = "36",
#'   mainCellType = "MC2",
#'   interactingCellType = "SC7",
#'   depedentMarker = "Podoplanin",
#'   interactive = FALSE,
#'   plotModelFit = FALSE,
#'   method = "lm",
#'   modelType = "dist200_"
#' )}
#'
#' @export
#' @rdname visualiseImageRelationship
#' @importFrom dplyr filter
#' @importFrom ggplot2
#'   ggplot scale_fill_distiller stat_density_2d geom_point theme_classic
#'   aes_string ggtitle facet_wrap aes xlab ylab ggtitle autoplot
#' @importFrom plotly ggplotly
#' @importFrom magrittr %>%
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
visualiseImageRelationship <- function(data,
                                       Rs,
                                       imageID,
                                       mainCellType,
                                       interactingCellType,
                                       depedentMarker,
                                       sizeVariable = NULL,
                                       shape = NULL,
                                       interactive = FALSE,
                                       plotModelFit = FALSE,
                                       method = "lm",
                                       modelType = "dist200_") {
  
  metadata_name <- paste0("Rs", Rs, "")
  
  SCE <- data 
  data <- metadata(SCE)[[metadata_name]]
  
  if(!depedentMarker %in% colnames(data)) {
    stop("The argument depedentMarker needs to exist in the data")
  }
  if(!imageID %in% data$imageID) {
    stop("The argument imageID needs to exist in the data")
  }
  
  data <- data[data$imageID == imageID, ]
  data$OriginalMarker <- data[, depedentMarker, drop = TRUE]
  data$fittedValues <- NA
  
  relationshipFormula <- paste0(
    depedentMarker, "~", paste0(modelType, interactingCellType) 
  )
  modelData <- data[data$cellType == mainCellType, ]
  model <- lm(formula(relationshipFormula), modelData)
  
  print(summary(model))
  
  data[data$cellType == mainCellType, "fittedValues"] <- predict(
    model,
    data[data$cellType == mainCellType, ]
  )
  
  
  if (plotModelFit == TRUE) {
    data[data$cellType == mainCellType, depedentMarker] <- predict(
      model, data[data$cellType == mainCellType, ]
    )
  }
  
  g1 <- ggplot2::ggplot() +
    ggplot2::stat_density_2d(
      data = data[data$cellType == interactingCellType, ],
      ggplot2::aes(
        x = x, y = y, fill = ..density..
      ),
      geom = "raster",
      contour = FALSE
    ) +
    ggplot2::scale_fill_distiller(palette = "Reds", direction = 1) +
    ggplot2::geom_point(
      data = data[data$cellType == mainCellType, ],
      ggplot2::aes_string(
        x = "x", y = "y",
        colour = depedentMarker,
        size = sizeVariable,
        shape = shape
      )
    ) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(paste(
      "Cell Points:", mainCellType, ",",
      "Cell Density:", interactingCellType, ",",
      "Model Fit:", plotModelFit
    )) +
    ggplot2::facet_wrap(~imageID, scales = "free")
  
  
  g2 <- data %>%
    dplyr::filter(cellType == mainCellType) %>%
    ggplot2::ggplot(
      ggplot2::aes_string(
        x = paste0(modelType, interactingCellType),
        y = "OriginalMarker"
      )
    ) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = lm) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("State Change Scatter Plot") +
    ggplot2::ylab(depedentMarker) +
    ggplot2::ylim(-1, NA)
  
  
  
  g3 <- data %>%
    dplyr::filter(cellType == mainCellType) %>%
    ggplot2::ggplot(ggplot2::aes_string(
      x = "OriginalMarker",
      y = "fittedValues"
    )) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = lm) +
    ggplot2::theme_classic() +
    ggplot2::xlab("True Values") +
    ggplot2::ylab("Fitted Values") +
    ggplot2::ggtitle("Predicted vs Real Values")
  
  # g4 <- ggplot2::autoplot(model) + ggplot2::theme_classic()
  
  if (interactive == TRUE) {
    g1 <- plotly::ggplotly(g1)
    g2 <- plotly::ggplotly(g2)
    g3 <- plotly::ggplotly(g3)
  }
  list(g1, g2, g3)
  # list(g1, g2, g3, g4)
}


