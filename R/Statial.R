#' Calculate Pairwise Distance Between Cell Types
#'
#' Calculates distance from each cell to the nearest cell of every type
#'
#' @param singleCellData
#' @param maxRS
#'
#' @examples
#' XYZ
#' 
#' @export
#' @rdname distanceCalculator
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select mutate inner_join arrange group_by slice rename relocate full_join
#' @importFrom spatstat.geom owin ppp closepairs
#' @importFrom tidyr pivot_wider
distanceCalculator = function(singleCellData, maxRS = 200){
  
  singleCellData = singleCellData %>% 
    tibble::rownames_to_column("i") %>% 
    dplyr::select(-i) %>% 
    tibble::rownames_to_column("cellIndex") %>% 
    dplyr::mutate(cellIndex = as.numeric(cellIndex)) 
  
  ow = spatstat.geom::owin(xrange = range(singleCellData$x), yrange = range(singleCellData$y))
  pppData = spatstat.geom::ppp(
    x = singleCellData$x ,
    y = singleCellData$y,
    window = ow,
    marks = singleCellData$cellType
  )
  
  closePairData = spatstat.geom::closepairs(pppData, rmax = maxRS)
  distanceData = data.frame(cellIndexA = closePairData$i, cellIndexB = closePairData$j, d = closePairData$d)
  
  
  cellAInformation = singleCellData %>% 
    dplyr::select(cellIndexA = cellIndex, cellTypeA = cellType, cellIDA = cellID) 
  cellBInformation = singleCellData %>% 
    dplyr::select(cellIndexB = cellIndex, cellTypeB = cellType) 
  
  processedDistanceData = distanceData %>% 
    dplyr::inner_join(cellAInformation, by = "cellIndexA") %>% 
    dplyr::inner_join(cellBInformation,by = "cellIndexB") %>% 
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
  
  
  processedDistanceData = processedDistanceData %>% 
    # XYZ - Not sure if this should be a full join
    dplyr::full_join(singleCellData, by = c("imageID", "cellID", "cellType"))
  
}


#' Calculates Pairwise Distance Between Cell Types Per Image
#'
#' Calculates distance from each cell to the nearest cell of every type per image
#'
#' @param singleCellData
#' @param nCores
#' @param Rs
#'
#' @examples
#' XYZ
#' 
#' @export
#' @rdname calculateCellDistances
#' @importFrom dplyr bind_rows mutate_if mutate_at rename_with contains starts_with vars
#' @importFrom spatstat.geom owin ppp closepairs
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom purrr reduce
#' @importFrom stringr str_replace
calculateCellDistances = function(singleCellData, splitBy = "imageID", nCores = 1, Rs = c(200)){
  singleCellDataDistances = singleCellData %>%
    split(~ imageID) %>% 
    BiocParallel::bplapply(distanceCalculator, rmax = max(Rs), BPPARAM  = BiocParallel::MulticoreParam(workers = nCores)) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate_if(is.numeric,function(x) ifelse(is.infinite(x), NA, x)) %>% 
    lapply(Rs, function(x) x %>% 
                             dplyr::mutate_at(dplyr::vars(dplyr::contains('dist_')), function(x) ifelse(x <= rmax, x, NA )) %>% 
                             dplyr::rename_with(function(x) stringr::str_replace(x, "dist_", paste0("dist", rmax, "_")), dplyr::starts_with("dist_"))
      , x = .) %>% 
    purrr::reduce(full_join)
}



#' Calculates Imhomogenous K Function Between Cell Types
#'
#' Calculates Imhomogenous K Function Between Cell Types by Image
#'
#' @param singleCellData
#' @param maxRS
#'
#' @examples
#' XYZ
#' 
#' @export
#' @rdname calculateK
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate_if bind_rows inner_join relocate 
#' @importFrom stringr word
#' @importFrom BiocParallel bplapply MulticoreParam
calculateK = function(singleCellData, nCores = 1, Rs = c(5,10, 20, 50, 100, 200, 400)){
  singleCellDataK = singleCellData %>%
    split(~ imageID) %>%
    BiocParallel::bplapply(lisaClust:::inhomLocalK,Rs = Rs, BPPARAM = BiocParallel::MulticoreParam(workers = nCores)) %>%
    lapply(as.data.frame) 
  
  # Correcting column names when Rs greater than max Rs of image
  correctedRadius = singleCellDataK %>% 
    lapply(function(x) unlist(lapply(colnames(x), function(x) min(Rs[Rs >= as.numeric(stringr::word(x, 1 ,sep = "_"))]))))
  cellTypeName = singleCellDataK %>% 
    lapply(function(x) word(colnames(x), 2, -1 ,sep = "_"))
  correctedColumnNames = correctedRadius %>% 
    mapply(function(newRadius, cellTypeName) mapply(function(a,b) paste0(a, "_",b), newRadius, cellTypeName), 
           newRadius = .,
           cellTypeName = cellTypeName,
           SIMPLIFY = FALSE)
  
  singleCellDataK = singleCellDataK %>% 
    mapply(function(x, newNames) x %>% purrr::set_names(newNames), x = ., newNames = correctedColumnNames, SIMPLIFY = FALSE) %>% 
    mapply(function(x, imageID) x %>%
             dplyr::rename_all( ~ paste0("count", .x)) %>%
             dplyr::mutate(imageID = imageID), x =., imageID = names(.), SIMPLIFY = FALSE) %>%
    dplyr::bind_rows() %>%
    tibble::rownames_to_column("cellID") %>%
    dplyr::mutate_if(is.numeric,function(x) ifelse(is.na(x), 0, x)) %>%
    # XYZ - Not sure if this should be a inner join
    dplyr::inner_join(singleCellData, by = c("imageID", "cellID")) %>% 
    dplyr::relocate(imageID)
}


#' Calculates Contamination Scores Based Using Random Forest
#'
#' Calculates contamination scores using random forest 
#'
#' @param singleCellData
#' @param  nCores 
#' @param seed
#' @param num.trees
#' @param verbose
#'
#' @examples
#' XYZ
#' 
#' @export
#' @rdname randomForestContaminationCalculator
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate select bind_rows inner_join  
#' @importFrom ranger ranger
#' @importFrom tibble column_to_rownames rownames_to_column
randomForestContaminationCalculator = function(singleCellData, 
                                               seed = 2022, 
                                               num.trees = 100,
                                               verbose = FALSE){
  rfData = singleCellData %>% 
    dplyr::mutate(cellID = paste0("cellID", cellID)) %>% 
    tibble::column_to_rownames("cellID") %>% 
    dplyr::select(cellType, markers) 
  
  set.seed(seed)
  rfModel <- ranger::ranger(as.factor(cellType) ~ ., data = rfData, num.trees = num.trees, probability = TRUE)
  
  if (verbose == TRUE){
    print(rfModel)
  }
  
  predictions = predict(rfModel, rfData)$predictions
  
  # XYZ - Put in Utility
  maxn <- function(n) function(x) order(x, decreasing = TRUE)[!is.na(x)][n]
  
  rfData = cbind(rfData, predictions) %>% 
    dplyr::mutate(rfMaxCellProb = apply(.[colnames(predictions)], 1, function(x) x[maxn(1)(x)])) %>% 
    dplyr::mutate(rfSecondLargestCellProb = apply(.[colnames(predictions)], 1, function(x) x[maxn(2)(x)])) %>% 
    dplyr::mutate(rfMainCellProb = apply(.[c("cellType", colnames(predictions))], 1, function(x) as.numeric(x[x["cellType"]]))) %>% 
    dplyr::select(-colnames(predictions)) %>% 
    tibble::rownames_to_column("cellID") %>% 
    dplyr::mutate(cellID = str_replace(cellID, "cellID", "")) 
  
  singleCellData = singleCellData %>% 
    dplyr::inner_join(rfData)
  
  singleCellData
  
  
}










#' Calculates Linear Models For Pairwise Interactions for all Variations of the Type Argument
#'
#' Calculates linear models For pairwise interactions for all variations of the type argument.
#'
#' @param singleCellData
#' @param markers
#' @param typeAll prefix to build models on
#' @param covariates
#' @param method 
#' @param isMixed
#' @param verbose
#' @param timeout
#' @param  nCores 
#'
#' @examples
#' XYZ
#' 
#' @export
#' @rdname calculateModelsAllInteractions
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom dplyr arrange group_by  summarise_at mutate_if mutate bind_rows left_join  
#' @importFrom tidyr gather
#' @importFrom janitor %>% 
calculateModelsAllInteractions = function(singleCellData, 
                                     markers, 
                                     typeAll = c("dist"), 
                                     covariates = NULL, 
                                     method = "lm", 
                                     isMixed = TRUE, 
                                     verbose = TRUE,
                                     timeout = 10,
                                     nCores = 1){
  
  typeVector = singleCellData %>% 
    dplyr::select(contains(typeAll)) %>% 
    colnames(.) %>% 
    stringr::str_split("_") %>% 
    lapply(function(x) x[[1]]) %>% 
    unlist() %>% 
    unique() %>% 
    paste0("_")

  allModels = typeVector %>% 
    mapply(calculateInteractionModels, type = ., MoreArgs = list(singleCellData = singleCellData,
                                                    markers = markers,
                                                    covariates = covariates, 
                                                    nCores = nCores,
                                                    method = method,
                                                    isMixed = isMixed, 
                                                    verbose = verbose,
                                                    timeout = timeout,
                                                    modelType), SIMPLIFY = FALSE) %>% 
    dplyr::bind_rows()
  
}




#' Calculates Linear Models For Pairwise Interactions
#'
#' Calculates linear models For pairwise interactions with options of mixed models and robust regression. 
#'
#' @param singleCellData
#' @param markers
#' @param type prefix to build models on
#' @param covariates
#' @param method 
#' @param isMixed
#' @param verbose
#' @param timeout
#' @param  nCores 
#'
#' @examples
#' XYZ
#' 
#' @export
#' @rdname calculateInteractionModels
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom dplyr arrange group_by  summarise_at mutate_if mutate bind_rows left_join  
#' @importFrom tidyr gather
#' @importFrom janitor %>% 
calculateInteractionModels = function(singleCellData, 
                                markers, 
                                type = "dist200_", 
                                covariates = NULL, 
                                method = "lm",
                                isMixed = FALSE, 
                                verbose = TRUE,
                                timeout = 10,
                                nCores = 1){
  


  if (isMixed == TRUE){
    splitData = split(singleCellData, ~cellType)
  }else{
    splitData = split(singleCellData, ~imageID)
  }
  
  singleCellDataCellInteractionModels = splitData %>% 
    BiocParallel::bplapply(buildModelsByCellType,
                           markers = markers,
                           covariates = covariates,
                           type = type,
                           method = method,
                           isMixed = isMixed,
                           verbose = verbose,
                           timeout = timeout,
                           BPPARAM = BiocParallel::MulticoreParam(workers = nCores)) %>% 
    bind_rows() 


  singleCellDataCellInteractionModels = singleCellDataCellInteractionModels %>%
    dplyr::arrange(pValue) %>% 
    mutate(type = type) %>% 
    mutate(covariates = covariates) %>% 
    mutate(covariates = ifelse(is.null(covariates), "None", covariates))
  
  relativeExpressionData = singleCellData %>% 
    dplyr::group_by(independent = cellType) %>% 
    dplyr::summarise_at(markers, mean, na.rm = TRUE) %>% 
    dplyr::mutate_if(is.numeric, heatmaply::normalize) %>% 
    tidyr::gather(-independent, key = "dependent", value = "relativeExpression") %>% 
    dplyr::mutate(interactingCell = independent) %>% 
    dplyr::mutate(independent = paste0(type, independent))
  
  singleCellDataCellInteractionModels = singleCellDataCellInteractionModels %>% 
    dplyr::mutate(interactingCell = str_replace(independent, type, "")) %>% 
    dplyr::left_join(relativeExpressionData)
  
  singleCellDataCellInteractionModels
}


#' Initiate Model Building By Cell Type
#'
#' Initiate Model Building By Cell Type.
#'
#' @param x
#' @param markers
#' @param type prefix to build models on
#' @param covariates
#' @param method 
#' @param isMixed
#' @param verbose
#' @param timeout
#' @param  nCores 
#'
#' @examples
#' XYZ
#' 
#' @export
#' @rdname buildModelsByCellType
#' @importFrom dplyr group_by  group_modify ungroup 
#' @importFrom janitor %>% 
buildModelsByCellType = function(x, markers, type, covariates, method, isMixed, verbose, timeout){
  finalModelsOutputs = x %>% 
    dplyr::group_by(cellType) %>% 
    dplyr::group_modify(~ .x %>% modelsPerCell(markers, type, covariates, method, isMixed, verbose, timeout)) %>% 
    dplyr::ungroup()
}




#' Formula Builder and Model Runner
#'
#' Formula Builder and Model Runner
#'
#' @param x
#' @param markers
#' @param type prefix to build models on
#' @param covariates
#' @param method 
#' @param isMixed
#' @param verbose
#' @param timeout
#'
#' @examples
#' XYZ
#' 
#' @export
#' @rdname modelsPerCell
#' @importFrom dplyr bind_rows   
#' @importFrom stringr str_detect str_replace
#' @importFrom janitor %>% 
modelsPerCell = function(x, markers, type, covariates, method, isMixed, verbose, timeout){
  
  if (verbose == TRUE){
    print(unique(x$imageID))
    print(unique(x$cellType))
  }
  
  cells = colnames(x)[Reduce("|", mapply(function(x, type) stringr::str_detect(colnames(x), type), type, MoreArgs = list(x = x), SIMPLIFY = FALSE))] 
  
  # Univariate Models
  if (length(type) == 1){
    formulas = lapply(markers, function(x, c) paste(x, "~", c), c = cells) %>% 
      unlist()
    
  # Interaction Models - two types
  }else{
    for (i in type){
      cells = stringr::str_replace(cells, i, "")
    }
    cells = unique(cells)
    cells = lapply(cells, function(x) unlist(lapply(type, paste0, x))) %>% 
      lapply(function(x) paste0(c(paste0(x, collapse = " + "), paste0(x, collapse = ":")), collapse = " + ")) %>%
      unlist()
    formulas = lapply(markers, function(x, c) paste(x, "~", c), c = cells) %>% 
      unlist()
  }
  
  modelOutputs = formulas %>% 
    mapply(fitInteractionModels, 
           f = ., 
           MoreArgs = list(x = x,
                           covariates = covariates, 
                           method = method, 
                           isMixed = isMixed, 
                           timeout = timeout), 
           SIMPLIFY = FALSE) %>% 
    dplyr::bind_rows() 
}





#' Formula Builder and Model Runner
#'
#' Formula Builder and Model Runner
#'
#' @param x
#' @param markers
#' @param type prefix to build models on
#' @param covariates
#' @param method 
#' @param isMixed
#' @param verbose
#' @param timeout
#'
#' @examples
#' XYZ
#' 
#' @export
#' @rdname fitInteractionModels
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
#' @importFrom janitor %>% 
fitInteractionModels = function(x, f, covariates, method, isMixed, timeout = timeout){
  dependent = stringr::str_split(f, " ~ ")[[1]][1]
  independent = stringr::str_split(f, " ~ ")[[1]][2]
  
  independentCheck = stringr::str_split(independent, "\\:") %>%
    unlist() %>% 
    stringr::str_split(" \\+ ") %>% 
    unlist() %>% 
    unique()
  
  f = paste(c(f, covariates), collapse = " + ")
  independentSplit = unlist(stringr::str_split(independent, " \\+ "))
  
  # XYZ - Should I have timeout for everything?
  outputs = try({
    
    if (isMixed == TRUE){
      if (modelType == "rlm"){
        
        model = R.utils::withTimeout({robustlmm::rlmer(formula(paste(f, "+ (1|patient/imageID)")), 
                                                       method = "DASvar",
                                                       data = x)},
                                     timeout = timeout, onTimeout = "silent")
        coefs = coef(summary(model))
        coefs = coefs[,colnames(coefs) != "df"]
        coefs = cbind(coefs, parameters::p_value(model)$p)
        modelSummary = summary(model)
        modelSummary$r.squared = NA
      }else{
        model = R.utils::withTimeout({lmerTest::lmer(formula(paste(f, "+ (1|patient/imageID)")), data = x)},
                                     timeout = timeout, onTimeout = "silent")
        coefs = coef(summary(model))
        coefs = coefs[,colnames(coefs) != "df"]
        modelSummary = summary(model)
        modelSummary$r.squared = MuMIn::r.squaredGLMM(model)[1]
      }
      
    }else{
      
      if (method == "rlm"){
        model = R.utils::withTimeout({MASS::rlm(formula(f), data = x)}, timeout = timeout, onTimeout = "silent")
        coefs = coef(summary(model))
        pValues = unlist(lapply(independentSplit, function(x) sfsmisc::f.robftest(model, var = x)$p.value))
        names(pValues) = independentSplit
        coefs = cbind(coefs, data.frame(pValue = NA))
        coefs[names(pValues), "pValue"] = pValues
        modelSummary = summary(model)
      }else{
        model = lm(formula(f), data = x)
        coefs = coef(summary(model))
        modelSummary = summary(model)
      }
    }
    
    beta = coefs[independentSplit,1]
    tValue =  coefs[independentSplit,3]
    pValue  =  coefs[independentSplit,4]
    rValue = modelSummary$r.squared
    sampleSize = length(modelSummary$residuals)
    outputs = data.frame(beta = beta, 
                         tValue = tValue,
                         pValue = pValue,
                         independent = independentSplit, 
                         dependent = dependent,
                         rValue =  rValue, 
                         sampleSize = sampleSize,
                         formula = f, 
                         isSingular = performance::check_singularity(model)
                         ) 
    
    if (isMixed == FALSE){
      outputs = outputs %>% 
        mutate(imageID = unique(x$imageID))
    }
    
  }, silent = TRUE)
  
  if (any(class(outputs) == "try-error")){
    outputs = data.frame(independent = independentSplit, 
                         dependent = dependent,
                         formula = f) 
    
    if (isMixed == FALSE){
      outputs = outputs %>% 
        mutate(imageID = unique(x$imageID))
    }
    
  }
  
  outputs
  
}



#' Convert Image Model Output to Cross-Validation Format
#'
#' Convert image model output to cross-validation format.
#'
#' @param imageModels
#' @param values_from
#' @param removeColsThresh
#' @param replacementValue 
#'
#' @examples
#' XYZ
#' 
#' @export
#' @rdname imageModelsCVFormat
#' @importFrom dplyr bind_rows mutate_if mutate rename select
#' @importFrom stringr str_detect str_replace str_split
#' @importFrom tidyr pivot_wider
#' @importFrom janitor %>% 
imageModelsCVFormat = function(imageModels, values_from = "tValue",  removeColsThresh = 0.2, replacementValue = 0){
  cvData = imageModels %>% 
    dplyr::mutate_if(is.numeric, function(x) ifelse(is.finite(x), x, NA)) %>% 
    dplyr::mutate(relationship = paste0(cellType, "_", dependent, "_", independent)) %>% 
    dplyr::rename(values_from = values_from) %>% 
    dplyr::select(imageID, relationship, values_from) %>% 
    tidyr::pivot_wider(names_from = relationship, values_from = values_from)
  
  cvData = cvData[,colSums(is.na(cvData)) < nrow(cvData)*removeColsThresh]
  cvData[is.na(cvData)] = replacementValue
  
  cvData
  
}