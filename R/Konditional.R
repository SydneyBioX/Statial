#' Identify changes in cell state of a cell type as it becomes closer to another, conditional on a third population.
#'
#' @param sce A SingleCellExperiment or SpatialExperiment
#'
#' @return A Statial result object
#'
#' @examples
#' 
#' 1+1
#' 
#' @export Konditional
#' @rdname Konditional
#' @importFrom BiocParallel SerialParam bplapply MulticoreParam
Konditional <- function(sce){
  1+1
}


#' Core function used by Konditional
#'
#'
#' @param image
#' @param r
#' @param from
#' @param to
#' @param parent
#' @param inhom
#' @param edge
#' @param includeOriginal
#' @param ... Any arguments passed into InhomLParent
#'
#' @return A single data frame row containing the konditional and orignal L values
#'
#' @examples
#' XYZ
#' 
#' @export
#' @rdname Konditional
#' @import spatstat
KonditionalCore = function(image,
                           r,
                           from,
                           to,
                           parent,
                           inhom = TRUE,
                           edge = FALSE,
                           includeOriginal = TRUE,
                           ...) {
    
    # Returns NA if to and from cell types not in image
    if(!(c(to, from)  %in% image$marks %>% all())) {
        condL = data.frame(child = NA, parent = NA)
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
            parent = unique(image$marks),
            edge = edge,
            inhom = inhom,
            ...
        )
    
    #return data frame of original and konditional values.
    condL = data.frame(original = originalVal,
                       konditional = konditionalVal)
    
    return(condL)
}
