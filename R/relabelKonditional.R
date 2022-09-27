
#' Function which randomises specific cells in an image and cacluates the Konditional value.
#'
#'
#' @param image A single image from a Single Cell Experiment object. 
#' @param nSim Number of randomisations which will be calculated. 
#' @param r Radius to evaluated pairwise relationships between from and to cells.
#' @param from The first cell type to be evaluated in the pairwise relationship.
#' @param to The second cell type to be evaluated in the pairwise relationship.
#' @param parent The parent population of the from cell type (must include from cell type).
#' @param inhom A logical value indicating whether to perform an inhomogeneous L function.
#' @param edge A logical value indicating whether to perform edge correction.
#' @param cores Number of cores for parallel processing.
#' @param ... Any arguments passed into Konditional
#' @return A data frame containing all images 
#'
#' @examples
#`
#' 
#' 
#' @export
#' @rdname relabelKonditional
#' @import dplyr
#' @import BiocParallel
#' @import tidyverse

relabelKonditional = function(image,
                              nSim = 1,
                              r,
                              from,
                              to,
                              parent,
                              inhom = TRUE,
                              edge = FALSE,
                              cores = 1,
                              ...
) {
    
    imageArray = replicate(nSim, image, simplify = FALSE)
    
    
    #relabel cells in parent population for all images in imageArray
    relabeled = bplapply(imageArray,
                         relabel,
                         labels = parent,
                         BPPARAM = MulticoreParam(workers = cores))
    
    relabeled = c(list(image), relabeled)
    names(relabeled) = seq_len(length(relabeled))
    
    parentDf = data.frame(from = from,
                          to = to, 
                          parent = I(list(parent)))
    
    #Calculated the child and parent values for the relabeled images
    relabeledDf = Konditional(
        relabeled,
        r = r,
        parentDf = parentDf,
        inhom = inhom,
        edge = edge,
        cores = cores,
        ...
    )
    
    relabeledDf = relabeledDf %>% 
        select(imageID, original, konditional) %>% 
        mutate(images = relabeled,
               type = ifelse(imageID == 1, "original", "randomised"))
    
    return(relabeledDf)
}



#' Function to randomise all specified cells labels in an image
#'
#'
#' @param image A single image from a Single Cell Experiment object. 
#' @param labels A vector of CellTypes to be randomised. If NULL all cells labels will be radomised.
#'
#' @return A data frame containing all pairwise cell relationships and their corresponding parent
#'
#' @examples
#`
#' 
#' 
#' @export
#' @rdname relabelKonditional
#' @import dplyr
#' 
relabel = function(image, labels = NULL) {
    
    #if labels are NULL relabel the whole image, otherwise relabel just the specified marks
    if (is.null(labels)) {
        relabeledCells = image %>% mutate(cellType = sample(cellType))
    } else {
        #split up cells into subset which will be relabeled, and subset which wont be relabeled
        notRelabel = image %>% filter(!(cellType %in% labels))
        toRelabel = image %>% filter(cellType %in% labels)
        
        #relabel toRelabel cells
        relabeled = toRelabel %>% mutate(cellType = sample(cellType))
        
        #combine relabeled cells and non relabeled cells
        relabeledCells = rbind(relabeled, notRelabel)
    }
    
    return(relabeledCells)
}

