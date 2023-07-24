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
#' @param typeAll
#'   A prefix that appears on the column names of all cell state modelling
#'   covariates. The default value is "dist"
#' @param covariates
#'   A list of additional covariates to be included in the model being built.
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
#' intensitiesData <- data.frame(t(
#'   SummarizedExperiment::assay(headSCE, "intensities")
#' ))
#' spatialData <- data.frame(SummarizedExperiment::colData(headSCE))
#' markersToUse <- colnames(intensitiesData)
#' singleCellData <- cbind(
#'   spatialData[rownames(intensitiesData), ], intensitiesData
#' )
#' singleCellData <- singleCellData %>%
#'   mutate(
#'     across(all_of(markersToUse), function(x) ifelse(is.na(x), 0, x))
#'   ) %>%
#'   mutate(across(where(is.factor), as.character))
#'
#' singleCellDataDistances <- getDistances(singleCellData,
#'   nCores = 1,
#'   Rs = c(200),
#'   whichCellTypes = c("MC2", "SC7")
#' )
#'
#' imageModels <- testStateChanges(
#'   singleCellData = singleCellDataDistances,
#'   markers = markersToUse,
#'   typeAll = c("dist200"),
#'   cellTypesToModel = "MC2",
#'   nCores = 1
#' )
#'
#' @export
#' @rdname testStateChanges
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom dplyr
#'   arrange group_by  summarise_at mutate bind_rows left_join filter
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
testStateChanges <- function(singleCellData,
                            Rs,
                            typeAll = c("dist"),
                            covariates = NULL,
                            method = "lm",
                            isMixed = FALSE,
                            randomIntercepts = c("imageID"),
                            cellTypesToModel = NULL,
                            verbose = FALSE,
                            timeout = 10,
                            nCores = 1) {
  
  markers = rownames(singleCellData)
  
  metadata_name <- paste0("dist", Rs, "")
  
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
                                 method = "lm",
                                 isMixed = FALSE,
                                 randomIntercepts = c("imageID"),
                                 verbose = TRUE,
                                 timeout = 10,
                                 nCores = 1) {
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
                                  method,
                                  isMixed,
                                  randomIntercepts,
                                  verbose,
                                  timeout) {
  finalModelsOutputs <- subsettedSingleCellData %>%
    dplyr::group_by(cellType) %>%
    dplyr::group_modify(
      ~ .x %>% modelsPerCellType(
        markers, type, covariates, method, isMixed,
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
  independentSplit <- unlist(stringr::str_split(independent, " \\+ "))
  
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
      
      beta <- coefs[independentSplit, 1]
      tValue <- coefs[independentSplit, 3]
      pValue <- coefs[independentSplit, 4]
      rValue <- modelSummary$r.squared
      sampleSize <- length(modelSummary$residuals)
      outputs <- data.frame(
        beta = beta,
        tValue = tValue,
        pValue = pValue,
        independent = independentSplit,
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
      independent = independentSplit,
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