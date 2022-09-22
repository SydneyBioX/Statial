#' Creates simulated image of tumour immune compartmentalisation
#'
#'
#' @param removal
#' @param r
#' @param sigma
#' @param k
#' @param mu
#' @param childSigma
#' @param includeTissue
#'
#'
#' @return A list containing two images, insig: where there is no conditional relationship, and sig: where there is a conditional relationship
#'
#'
#' @examples
#' images = simulateCompartment()
#' plot(images$sig)
#' 
#' @export
#' @rdname simulateCompart
#' @import spatstat
#' @importFrom EBImage gblur
simulateCompartment = function (removal = 0.25,
                                r = 0.1,
                                sigma = 0.05,
                                k = 40,
                                mu = 50,
                                childSigma = 1,
                                includeTissue = TRUE) {
    
    #constructing compartment densities (cDen) which other densities will be based off
    compartment = rMatClust(kappa = k, r = r, mu = mu)
    cDen = density(compartment, sigma = sigma)
    
    #Defining tumour den
    tumourDen = cDen
    
    #Density values in the bottom removal% will = 0, so that no points can be placed there, make top 1-removal % be solid (create a mask)
    tumourDen[tumourDen < max(tumourDen) * removal] = 0
    tumourDen[tumourDen > 0] = max(tumourDen)
    
    #Invert tumourDen mask to create T cell mask, doing pmax with zero so that there are no negative probabilities 
    tDen = ((-tumourDen) / sum(tumourDen) * sum(tumourDen)) + max(tumourDen)
    tDen$v = pmax(tDen$v, 0)
    
    #defining Cd8 densities
    cd8Den = cDen
    cd8Den = (tDen * cDen) / mean(max(cDen), max(tDen))
    cd8Den = ((1 - cd8Den) / sum(cd8Den) * sum(cd8Den)) + max(cd8Den)
    cd8Den[cd8Den > max(cd8Den) * removal] = 0
    
    #Smooth out cd8, use pmax so there are no negative probabilities, and scale the values up so that the cell counts matches the other densities
    cd8Den$v = gblur(cd8Den, sigma = childSigma)
    cd8Den$v = pmax(cd8Den$v, 0)
    cd8Den = cd8Den * (mean(max(tumourDen), max(tDen)) / max(cd8Den))
    
    #Make cells using densities
    tumourCells = rpoispp(tumourDen)
    tCells = rpoispp(tDen)
    cd8Insig = rpoispp(tDen)
    cd8Sig = rpoispp(cd8Den)
    
    #Define marks
    marks(tumourCells) = factor("tumour_cells")
    marks(tCells) = factor("t_cells")
    marks(cd8Insig) = factor("cd8_t_cells")
    marks(cd8Sig) = factor("cd8_t_cells")
    
    #set up simulations.
    simInsig = superimpose(tumourCells, tCells, cd8Insig)
    simSig = superimpose(tumourCells, tCells, cd8Sig)
    
    if (includeTissue == TRUE) {
        tissueCells = rpoispp(mean(intensity(tumourCells), intensity(tCells)))
        marks(tissueCells) = factor("tissue_cells")
        simInsig = superimpose(simInsig, tissueCells)
        simSig = superimpose(simSig, tissueCells)
    }
    
    return(list("insig" = simInsig, "sig" = simSig))
}

