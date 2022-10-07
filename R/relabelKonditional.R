
#' Function which randomises specific cells in an image and cacluates the Konditional value.
#'
#'
#' @param image A single image from a Single Cell Experiment object. 
#' @param nSim Number of randomisations which will be calculated. 
#' @param r Radius to evaluated pairwise relationships between from and to cells.
#' @param from The first cell type to be evaluated in the pairwise relationship.
#' @param to The second cell type to be evaluated in the pairwise relationship.
#' @param parent The parent population of the from cell type (must include from cell type).
#' @param returnImages A logical value to indicate whether the function should 
#' return the randomised images along with the Konditional values. 
#' @param inhom A logical value indicating whether to account for inhomogeneity.
#' @param edge A logical value indicating whether to perform edge correction.
#' @param cores Number of cores for parallel processing.
#' @param ... Any arguments passed into \code{\link[Statial]{Konditional}}
#' @return A data frame containing Konditional value for each randomised image. 
#' If `returnImages = TRUE` function will return a list with Konditional values 
#' and the randomised images.
#'
#' @examples
#' data("exampleImage")
#' 
#' relabelResult = relabelKonditional(
#' image = exampleImage,
#' nSim = 5,
#' r = 0.05,
#' from = "cd8_t_cells",
#' to = "tumour_cells",
#' parent = c("cd8_t_cells", "t_cells"),
#' cores = 40)
#' 
#' @export
#' @rdname relabelKonditional
#' @import dplyr
#' @import BiocParallel
#' @import tidyverse

relabelKonditional <- function(cells,
                              nSim = 1,
                              r,
                              from,
                              to,
                              parent,
                              returnImages = FALSE,
                              inhom = TRUE,
                              edge = FALSE,
                              cores = 1,
                              ...
) {
    
    imageArray <- replicate(nSim, cells, simplify = FALSE)
    
    
    #relabel cells in parent population for all images in imageArray
    relabeled <- bplapply(imageArray,
                         relabel,
                         labels = parent,
                         BPPARAM = MulticoreParam(workers = cores))
    
    relabeled <- c(list(cells), relabeled)
    names(relabeled) <- as.character(seq_along(relabeled))
    
    parentDf <- data.frame(from = from,
                          to = to, 
                          parent = I(list(parent)))
    
    #Calculated the child and parent values for the relabeled images
    relabeledDf <- Konditional(
        imageData = relabeled,
        r = r,
        parentDf = parentDf,
        inhom = inhom,
        edgeCorrect = edge,
        cores = cores,
        ...
    )
    
    relabeledDf <- relabeledDf %>% 
        select(imageID, original, konditional, r) %>% 
        mutate(type = ifelse(imageID == 1, "original", "randomised"))
    
    if(returnImages) {
        return(list(relabeledDf, relabeled))
    }
    
    return(relabeledDf)
}



#' Function to permute all specified cells labels in an image
#'
#'
#' @param image A single image from a Single Cell Experiment object. 
#' @param labels A vector of CellTypes labels to be permuted If NULL all cells 
#' labels will be radomised.
#'
#' @return A data frame containing all pairwise cell relationships and their 
#' corresponding parent
#'
#' @examples
#' data("exampleImage")
#' 
#' #Permute CD8 T cells and T cell labels in the image
#' relabeledImage = relabel(exampleImage, labels = c("cd8_t_cells", "t_cells"))
#' plot(relabeledImage)
#' 
#' @export
#' @rdname relabelKonditional
#' @import dplyr
#' 
relabel <- function(image, labels = NULL) {
    
    #if labels are NULL relabel the whole image, otherwise relabel just the specified marks
    if (is.null(labels)) {
        relabeledCells <- image %>% mutate(cellType = sample(cellType))
    } else {
        #split up cells into subset which will be relabeled, and subset which wont be relabeled
        notRelabel <- image %>% filter(!(cellType %in% labels))
        toRelabel <- image %>% filter(cellType %in% labels)
        
        #relabel toRelabel cells
        relabeled <- toRelabel %>% mutate(cellType = sample(cellType))
        
        #combine relabeled cells and non relabeled cells
        relabeledCells <- rbind(relabeled, notRelabel)
    }
    
    return(relabeledCells)
}

