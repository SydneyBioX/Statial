#' Identify changes in cell state of a cell type as it becomes closer to another, conditional on a third population.
#'
#' @param imageData A SingleCellExperiment or SpatialExperiment or a list of single images. 
#' @param parentDf A data frame from \code{\link[Statial]{parentCombinations}}
#' @param r Radius to evaluated pairwise relationships between from and to cells.
#' @param from The first cell type to be evaluated in the pairwise relationship.
#' @param to The second cell type to be evaluated in the pairwise relationship.
#' @param parent The parent population of the from cell type (must include from cell type).
#' @param inhom  A logical value indicating whether to account for inhomogeneity.
#' @param edgeCorrect A logical value indicating whether to perform edge correction.
#' @param window Type of window for data, either `square`, `convex` or `concave`, passed into \code{\link[Statial]{makeWindow}}
#' @param window.length A tuning parameter for controlling the level of concavity when estimating concave windows.
#' Passed into \code{\link[Statial]{makeWindow}}
#' @param weightQuantile A decimal value indicating what quantile of parent density used to weight the `from` cells.
#' @param includeZeroCells A logical value indicating whether to include cells with
#' zero counts in the pairwise association calculation.
#' @param includeOriginal A logical value to return the original L function values along with the konditional values.
#' @param cores Number of cores for parallel processing.
#' @return A Koditional result object
#'
#' @examples
#' 1+1
#' 
#' @export Konditional
#' @rdname Konditional
#' @import dplyr
#' @import tidyr
#' @import BiocParallel

Konditional = function(imageData,
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
                       cores = 1)
{
    
    if(is.null(c(parentDf, from, to, parent))){
        stop("Please specificy a parentDf (obtained from parentCombinations, or from, to, and parent cellTypes")
    } else if(!is.null(from) &
              !is.null(to) &
              !is.null(parent)) {
        
        parentDf = data.frame(from = from,
                              to = to,
                              parent = I(list(parent)))
        
    }
    
    
    # Creating a vector for images
    if(class(imageData) == "SingleCellExperiment") {
        imageData = imageData %>%
            colData() %>%
            data.frame()
        
        imageData = mutate(imageData, imageID = as.character(imageID))
        images = split(imageData, imageData$imageID)
    }
    
    if(class(imageData) != "list") {
        imageData = list(imageData)
        names(imageData) = as.character(seq_along(imageData))
    }
    
    images = imageData
    
    imagesInfo = data.frame(imageID = names(images),
                            images = I(images))
    
    # Create all combinations of specified parameters
    allCombinations = expand_grid(
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
    konditionalDf = merge(imagesInfo, allCombinations, all = TRUE) %>% 
        mutate(test = paste(from, "__", to))
    
    if("parent_name" %in% names(konditionalDf)) {
        konditionalDf = mutate(konditionalDf, test = paste(test, "__", parent_name))
    }
    
    # Calculate conditional L values
    lVals = bpmapply(
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
        BPPARAM = MulticoreParam(workers = cores)
    )
    
    # Combine data.frame rows
    lVals = lVals %>%
        bind_rows()
    
    lValsClean = konditionalDf %>%
        mutate(parent_name = "") %>% 
        select(-c(images, from, to, parent_name, parent)) %>%
        cbind(lVals) %>%
        remove_rownames() 
    
    if(includeOriginal == FALSE) {
        lValsClean = lValsClean %>% 
            select(imageID, test, konditional, r, weightQuantile, inhom, edge, includeZeroCells, window.length)
    } else {
        lValsClean = lValsClean %>% 
            select(imageID, test, original, konditional, r, weightQuantile, inhom, edge, includeZeroCells, window.length) 
    }
    
    return(lValsClean)
}



#' Core function used by Konditional
#'
#'
#' @param image A single image from a Single Cell Experiment object. 
#' @param r The radius which pairwise cell relationships are evaluated at.
#' @param from The first cell type to be evaluated in the pairwise relationship.
#' @param to The second cell type to be evaluated in the pairwise relationship.
#' @param parent The parent population of the from cell type (must include from cell type).
#' @param inhom A logical value indicating whether to perform an inhomogeneous L function.
#' @param edge A logical value indicating whether to perform edge correction.
#' @param includeOriginal A logical value to return the original L function values along with the konditional values.
#' @param ... Any arguments passed into \code{\link[Statial]{inhomLParent}}.
#'
#' @return A single row of a data frame containing the konditional and orignal L values.
#'
#' @examples
#' XYZ
#' 
#' @export
#' @rdname KonditionalCore
#' @import spatstat
#' @import tidyverse
KonditionalCore = function(image,
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
    if(!(c(to, from)  %in% unique(image$cellType) %>% all())) {
        condL = data.frame(original = NA, konditional = NA)
        rownames(condL) = paste(from, "__", to)
        return(condL)
        
    }
    
    #Calculated the child and parent values for the image
    konditionalVal =
        inhomLParent(
            image,
            Rs = r,
            from = from,
            to = to,
            parent = parent,
            edge = edge,
            inhom = inhom,
            weightQuantile = weightQuantile,
            ...
        )
    
    
    if(includeOriginal == FALSE) {
        condL = data.frame(konditional = konditionalVal)
        return(condL)
    }
    
    originalVal =
        inhomLParent(
            image,
            Rs = r,
            from = from,
            to = to,
            parent = unique(image$cellType),
            edge = edge,
            inhom = inhom,
            weightQuantile = weightQuantile,
            ...
        )
    
    #return data frame of original and konditional values.
    condL = data.frame(original = originalVal,
                       konditional = konditionalVal)
    
    return(condL)
}


