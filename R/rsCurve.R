#' Evaluation of Konditional over a range of radii.
#'
#'
#' @param image A single image from a Single Cell Experiment object. 
#' @param from The first cell type to be evaluated in the pairwise relationship.
#' @param to The second cell type to be evaluated in the pairwise relationship.
#' @param parent The parent population of the from cell type (must include from cell type).
#' @param rs A vector of radii to evaluate konditional over.
#' @param inhom A logical value indicating whether to perform an inhomogeneous L function.
#' @param edge A logical value indicating whether to perform edge correction.
#' @param se A logical value to indicate if the standard deviation of konditional should be calculated to construct error bars.
#' @param nSim Number of randomisations to perform using \code{\link[Statial]{relabelKonditional}}, which will be used to calculated the SE.
#' @param cores Number of cores for parallel processing.
#' @param ... Any arguments passed into \code{\link[Statial]{inhomLParent}}.
#'
#' @return A data frame of original L values and Konditional values evaluated over a range of radii.
#'
#' @examples
#' XYZ
#' 
#' @export
#' @rdname rsCurve
#' @import tidyverse


rsCurve = function(image,
                   from,
                   to,
                   parent,
                   rs = seq(10, 100, 10),
                   inhom = TRUE,
                   edge = FALSE,
                   se = FALSE,
                   nSim = 20,
                   cores = 1,
                   ...) {
    
    
    konditionalVals = Konditional(imageData = image,
                                  from = from,
                                  to = to,
                                  parent = parent,
                                  r = rs,
                                  inhom = inhom,
                                  edge = edge,
                                  cores = cores,
                                  includeOriginal = TRUE,
                                  ...)
    
    rsDf = konditionalVals %>% 
        select(r, original, konditional)
    
    
    if (se == TRUE) {
        seDf = relabelKonditional(image = image,
                                  nSim = nSim,
                                  r = rs,
                                  from = from,
                                  to = to, 
                                  parent = parent,
                                  returnImages = FALSE,
                                  inhom = inhom,
                                  edge = edge,
                                  cores = cores,
                                  ...)
        
        seDf = seDf %>%
            filter(type != "original") %>%
            select(r, original, konditional) %>%
            group_by(r) %>%
            summarise(originalSd = sd(original),
                      konditionalSd = sd(konditional))
        
        rsDf = merge(rsDf, seDf, by = "r")
    }
    
    return(rsDf)
}