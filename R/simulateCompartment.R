#' Creates simulated image of tumour immune compartmentalisation
#'
#'
#' @param removal A decimal value indicating how much of the density map to 
#' discard to make a mask.
#' @param r Radius of the clusters passed to \code{\link[spatstat.random]{rMatClust}}.
#' @param sigma A numerical variable to indicated the standard deviation of the 
#' isotropic smoothing kernel for density calculation,
#' passed into \code{\link[spatstat.core]{density.ppp}}.
#' @param k Intensity of the cluster centers passed to \code{\link[spatstat.random]{rMatClust}}.
#' @param mu Average number of cells per cluster, passed to \code{\link[spatstat.random]{rMatClust}}
#' @param childSigma A numerical variable to indicated the standard deviation 
#' of the isotropic smoothing kernel for density calculation,
#' used to create the cd8 population of cells. Passed into \code{\link[EBImage]{gblur}}.
#' @param includeTissue A logical to include a uniformly distributed tissue cell type.
#'
#'
#' @return A list containing two images, insig: where there is no conditional 
#' relationship, and sig: where there is a conditional relationship
#'
#'
#' @examples
#' images = simulateCompartment()
#' simSig = images$sig
#' 
#' plot(x = simSig$x, y = simSig$y, col = simSig$cellType)
#' 
#' @export
#' @rdname simulateCompartment
#' @import spatstat
#' @importFrom spatstat.random rMatClust rpoispp
#' @importFrom spatstat.core density.ppp
#' @importFrom EBImage gblur filter2
simulateCompartment = function (removal = 0.25,
                                r = 0.1,
                                sigma = 0.05,
                                k = 40,
                                mu = 50,
                                childSigma = 1,
                                includeTissue = TRUE) {
    
    #constructing compartment densities (cDen) which other densities will be based off
    compartment <- spatstat.random::rMatClust(kappa = k, scale = r, mu = mu)
    cDen <- spatstat.core::density.ppp(compartment, sigma = sigma)
    
    #Defining tumour den
    tumourDen <- cDen
    
    #Density values in the bottom removal% will = 0, so that no points can be placed there, make top 1-removal % be solid (create a mask)
    tumourDen[tumourDen < max(tumourDen) * removal] = 0
    tumourDen[tumourDen > 0] = max(tumourDen)
    
    #Invert tumourDen mask to create T cell mask, doing pmax with zero so that there are no negative probabilities 
    tDen <- ((-tumourDen) / sum(tumourDen) * sum(tumourDen)) + max(tumourDen)
    tDen$v <- pmax(tDen$v, 0)
    
    #defining Cd8 densities
    cd8Den <- cDen
    cd8Den <- (tDen * cDen) / mean(max(cDen), max(tDen))
    cd8Den <- ((1 - cd8Den) / sum(cd8Den) * sum(cd8Den)) + max(cd8Den)
    cd8Den[cd8Den > max(cd8Den) * removal] = 0
    
    #Smooth out cd8, use pmax so there are no negative probabilities, and scale the values up so that the cell counts matches the other densities
    cd8Den$v <- EBImage::gblur(cd8Den$v, sigma = childSigma)
    cd8Den$v <- pmax(cd8Den$v, 0)
    cd8Den <- cd8Den * (mean(max(tumourDen), max(tDen)) / max(cd8Den))
    
    #Make cells using densities
    tumourCells <- spatstat.random::rpoispp(tumourDen)
    tCells <- spatstat.random::rpoispp(tDen)
    cd8Insig <- spatstat.random::rpoispp(tDen)
    cd8Sig <- spatstat.random::rpoispp(cd8Den)
    
    #Define marks
    marks(tumourCells) <- factor("tumour_cells")
    marks(tCells) <- factor("t_cells")
    marks(cd8Insig) <- factor("cd8_t_cells")
    marks(cd8Sig) <- factor("cd8_t_cells")
    
    #set up simulations.
    simInsig <- superimpose(tumourCells, tCells, cd8Insig)
    simSig <- superimpose(tumourCells, tCells, cd8Sig)
    
    if (includeTissue == TRUE) {
        tissueCells <- spatstat.random::rpoispp(mean(intensity(tumourCells),
                                                     intensity(tCells)))
        marks(tissueCells) <- factor("tissue_cells")
        simInsig <- superimpose(simInsig, tissueCells) 
        simSig <- superimpose(simSig, tissueCells)
    }
    
    simInsig <- simInsig %>% PPPdf()
    simSig <- simSig %>% PPPdf()
    
    return(list("insig" = simInsig, "sig" = simSig))
}