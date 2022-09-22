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

#' Create all combinations of cell type relationships from a list of parents
#'
#'
#' @param all A list of all cell types
#' @param ... Vectors of each parent population
#'
#' @return A data frame containing all pairwise cell relationships and their corresponding parent
#'
#' @examples
#' tcells = c("CD4", "CD8")
#' tissue = c("epithelial", "stromal")
#' allCells = c("tumour", tissue, tcells)
#' 
#' parentCombinations(all = allCells, tcells, tissue)
#' 
#' 
#' @export
#' @rdname Konditional
#' @import dplyr
#' @import tidyr
parentCombinations = function(all, ...) {
    
    #Gets variable names of all the parent vector
    names = deparse(substitute(c(...))) %>%
        str_remove_all("^[c(]{2}|[)]$") %>% 
        str_split(", ") %>% 
        unlist()
    
    parentList = list(...)
    names(parentList) = names
    
    #Creates data.frame of parent name and parent vector
    parentTable = data.frame(parent_name = names(parentList),
                             parent = I(parentList))
    
    
    #Creates all combination of parent and child
    parentDfs = list()
    
    for (i in 1:length(names(parentList))) {
        parentDfs[[i]] =  crossing(from = parentList[[i]], parent_name = names(parentList)[i])
    }
    
    parentDf = bind_rows(parentDfs) %>% 
        merge(parentTable, by = "parent_name") %>% 
        expand_grid(to = unique(all)) %>% 
        data.frame() %>% 
        select(from, to, parent, parent_name) %>% 
        filter(from != to)
    
    
    
    return(parentDf)
    
}

#############
# SIMULATION FILE

#############

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








##########
#MAKE WINDOW FILE

#########

#' Creates a window for PPP object
#'
#'
#' @param data
#' @param window
#' @param window.length
#'
#' @examples
#' XYZ
#' 
#' @export
#' @rdname makeWindow
#' @import spatstat.geom
#' @import concaveman

makeWindow <- function(data,
                       window = "square",
                       window.length = NULL) {
    data = data.frame(data)
    ow <-
        spatstat.geom::owin(xrange = range(data$x), yrange = range(data$y))
    
    if (window == "convex") {
        p <- spatstat.geom::ppp(data$x, data$y, ow)
        ow <- spatstat.geom::convexhull(p)
        
    }
    if (window == "concave") {
        message(
            "Concave windows are temperamental. Try choosing values of window.length > and < 1 if you have problems."
        )
        if (is.null(window.length) | is.na(window.length)) {
            window.length <- (max(data$x) - min(data$x)) / 20
        } else{
            window.length <- (max(data$x) - min(data$x)) / 20 * window.length
        }
        dist <- (max(data$x) - min(data$x)) / (length(data$x))
        bigDat <-
            do.call("rbind", lapply(as.list(as.data.frame(t(data[, c("x", "y")]))), function(x)
                cbind(
                    x[1] + c(0, 1, 0,-1,-1, 0, 1,-1, 1) * dist,
                    x[2] + c(0, 1, 1, 1,-1,-1,-1, 0, 0) * dist
                )))
        ch <-
            concaveman::concaveman(bigDat,
                                   length_threshold = window.length,
                                   concavity = 1)
        poly <- as.data.frame(ch[nrow(ch):1, ])
        colnames(poly) <- c("x", "y")
        ow <-
            spatstat.geom::owin(
                xrange = range(poly$x),
                yrange = range(poly$y),
                poly = poly
            )
        
    }
    ow
}


##########
#INHOM L PARENT FILE

#########


#' Calculates conditional L value
#'
#'
#' @param data
#' @param Rs
#' @param window
#' @param window.length
#' @param weightQuantile
#' @param from
#' @param to
#' @param edgeCorrect
#' @param includeZeroCells
#' @param parent
#' @param inhom
#' @param closePairs
#'
#' @examples
#' XYZ
#' 
#' @export
#' @rdname inhomLParent
#' @import spatstat
inhomLParent <- function (data,
                          Rs = 20,
                          window = "convex",
                          window.length = NULL,
                          weightQuantile = .80,
                          from = NULL,
                          to = NULL,
                          edgeCorrect = TRUE,
                          includeZeroCells = TRUE,
                          parent = NULL,
                          inhom = TRUE,
                          closePairs = NULL) {
    
    if (class(data) == "ppp") {
        data = PPPdf(data)
    }
    
    #if class is data frame it show make window etc.
    ow <- makeWindow(data, window, window.length)
    
    #sigma = radius
    if (is.null(Rs))
        Rs = 50
    sigma = Rs
    
    X <-
        spatstat.geom::ppp(
            x = data$x,
            y = data$y,
            window = ow,
            marks = data$cellType
        )
    

    Y <- X
    if(!is.null(parent)) Y <- Y[Y$marks %in% parent, ]
    
    Area = area(X)
    
    
    if(!is.null(Rs)){
        den <- spatstat.core::density.ppp(Y, sigma = sigma, kernel = "disc") #Use disc kernel for interpretation
        den <- den / max(den)
        #den$v <- pmax(den$v, minLambda)
        Area = area(X)*mean(den) #area of the parent 
    }
    
    
    maxR <- min(ow$xrange[2]- ow$xrange[1], ow$yrange[2]- ow$yrange[1])/2.01
    Rs <- unique(pmin(c(0, sort(Rs)),maxR))
    
    
    if(is.null(from)) from <- levels(data$cellType)
    if(is.null(to)) to <- levels(data$cellType)
    
    use <- data$cellType %in% c(from, to)
    fulldata = data
    data <- data[use,]
    X <- X[use,]
    
    #pairwise relationships for a r
    if(is.null(closePairs)) {
        closePairs = spatstat.geom::closepairs(X, max(Rs), what = "ijd", distinct = FALSE) 
        closePairs$j <- data$cellID[closePairs$j]
        closePairs$i <- data$cellID[closePairs$i]
        
    } else {
        closePairs = closePairs %>% data.frame() %>% filter(d <= max(Rs))
        
        #Convert True false into indexes using seq_along and then subset closePairs dataframe
        
        closePairs = closePairs[closePairs$i %in% seq_along(use)[use],]
        closePairs = closePairs[closePairs$j %in% seq_along(use)[use],]
        
    }
    
    
    
    p <- closePairs
    
    n <- X$n
    
    cT <- data$cellType
    names(cT) <- data$cellID
    
    p$d <- cut(p$d, Rs, labels = Rs[-1], include.lowest = TRUE)
    
    w <- rep(1, length(X))
    p$wt <- rep(1,length(p$d))
    if(!is.null(sigma)){
        np <- spatstat.geom::nearest.valid.pixel(X$x, X$y, den)
        w <- den$v[cbind(np$row, np$col)]
        names(w) <- data$cellID
        invWeight = 1/w[as.character(p$i)]
        p$wt <- pmin(invWeight, quantile(invWeight, weightQuantile))
    }        
    
    num <- tapply(w, data$cellType, length)
    #weights of each cell type
    lam <- tapply(w, data$cellType, sum)/Area
    #sum of the weights 
    if(inhom)num <- tapply(pmin(1/w, quantile(1/w, weightQuantile)), data$cellType, sum)
    if(!inhom)p$wt <- 1
    
    # # inhom density
    # p$wt <- rep(1,length(p$d))
    # if(!is.null(sigma)){
    #     np <- spatstat.geom::nearest.valid.pixel(X$x, X$y, den)
    #     w <- den$v[cbind(np$row, np$col)]
    #     names(w) <- data$cellID
    #     p$wt <- w[p$j]#*mean(w)
    #     rm(np)
    # }
    
    
    #lam <- table(data$cellType)/spatstat.geom::area(X)
    # p$wt <- as.numeric(p$wt/lam[cT[p$j]])
    
    p$cellTypeJ <- cT[as.character(p$j)]
    p$cellTypeI <- cT[as.character(p$i)]
    p$i <- factor(p$i, levels = data$cellID)
    
    if(edgeCorrect){
        edge <- sapply(Rs[-1],function(x)borderEdge(X,x), simplify = FALSE)
        edge <- do.call("cbind",edge)
        edge <- as.data.frame(edge)
        colnames(edge) <- Rs[-1]
        edge$i <- data$cellID
        #edge$i <- factor(data$cellID, levels = data$cellID)
        edge <- tidyr::pivot_longer(edge,-i,"d")
        p <- dplyr::left_join(as.data.frame(p), edge, c("i", "d"))
    }else{
        p <- as.data.frame(p)
        p$value <- 1
    }
    
    
    p$d <- factor(p$d, levels = Rs[-1])
    
    p <- p[as.character(p$i) != as.character(p$j), ]
    
    use <- p$cellTypeI %in% from & p$cellTypeJ %in% to
    p <- p[use,]
    
    r <- inhomL(p, lam, X, Rs, num, Area)
    
    wt <- r$wt
    names(wt) <- paste(r$cellTypeI, r$cellTypeJ, sep = "__")
    
    m1 <- rep(from, times = length(to))
    m2 <- rep(to, each = length(from))
    labels <- paste(m1, m2, sep = "__")
    
    assoc <- rep(-sum(Rs), length(labels))
    names(assoc) <- labels
    if(!includeZeroCells)assoc[!(m1%in%X$marks&m2%in%X$marks)] = NA
    assoc[names(wt)] <- wt
    names(assoc) <- labels
    
    assoc
}





#' Calculates L value from weights and lambda values
#'
#'
#' @param p
#' @param lam
#' @param X
#' @param Rs
#' @param num
#' @param Area
#'
#' @examples
#' XYZ
#' 
#'
#' @rdname inhomLParent
#' @import data.table

inhomL <-
    function (p, lam, X, Rs, num, Area) {
        r <- data.table::as.data.table(p)
        
        r$wt <- (r$wt)/as.numeric(lam[r$cellTypeJ])/(as.numeric(num[r$cellTypeI]))

        r <- r[,j:=NULL]
        r <- r[,value:=NULL]
        r <- r[,i:=NULL]
        data.table::setkey(r, d, cellTypeI, cellTypeJ)
        r <- r[data.table::CJ(d, cellTypeI, cellTypeJ, unique = TRUE)
        ][, lapply(.SD, sum), by = .(d, cellTypeI, cellTypeJ)
        ][is.na(wt), wt := 0]
        r <- r[, wt := cumsum(wt), by = list(cellTypeI, cellTypeJ)]
        r <- r[, list(wt=sum(sqrt(wt/pi))), by=.(cellTypeI, cellTypeJ)]
        r$wt <- r$wt - sum(Rs)
        
        r <- as.data.frame(r)
        
        r
        
    }

#' Converts a PPP object to a data frame for inhomLParent function
#'
#'
#' @param ppp a PPP object
#'
#' @examples
#' XYZ
#' 
#' 
#' @rdname inhomLParent
#' @import spatstat.geom
#' 
PPPdf = function(ppp) {
    x <- as.data.frame(ppp)
    x$cellType = factor(x$marks)
    x$x <- x$x
    x$y <- x$y
    x$cellID <- factor(seq_len(nrow(x)))
    return(x)
}


#' Edge correction for L values
#'
#'
#' @param X
#' @param maxD
#'
#' @examples
#' XYZ
#' 
#' 
#' @rdname inhomLParent
#' @import spatstat.geom

borderEdge <- function(X, maxD){
    W <-X$window
    bW <- spatstat.geom::union.owin(spatstat.geom::border(W,maxD, outside = FALSE),
                                    spatstat.geom::border(W,2, outside = TRUE))
    inB <- spatstat.geom::inside.owin(X$x, X$y, bW)
    e <- rep(1, X$n)
    if(any(inB)){
        circs <-spatstat.geom:: discs(X[inB], maxD, separate = TRUE)
        circs <- spatstat.geom::solapply(circs, spatstat.geom::intersect.owin, X$window)
        areas <- unlist(lapply(circs, spatstat.geom::area))/(pi*maxD^2)
        e[inB] <- areas
    }
    
    e
}