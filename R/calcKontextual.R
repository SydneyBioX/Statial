#' @noRd
#' @importFrom methods is
#' @importFrom stats quantile

calcKontextual <- function(data,
                          child1,
                          child2,
                          parent,
                          r,
                          closePairs = NULL,
                          window = "convex",
                          window.length = NULL,
                          inhom = FALSE,
                          returnWeight = FALSE){
  
  if (is(data, "ppp")) {
    data <- PPPdf(data)
  }
  
  
  ####### !!!!!!! I dumped this because of errors  
  # if (!("cellID" %in% names(data))) {
  data$cellID <- factor(seq_len(nrow(data)))
  #  }
  
  # if class is data frame it show make window etc.
  ow <- Statial::makeWindow(data, window, window.length)
  
  if (is.null(r)) {
    r <- 50
  }
  
  # Filter data, !!!!!! wouldn't do this if looping later.  
  #  data <- data[data$cellType%in%c(child1, child2, parent),]
  
  # Create ppp  
  X <-
    spatstat.geom::ppp(
      x = data$x,
      y = data$y,
      window = ow,
      marks = data$cellType
    )
  
  Area <- spatstat.geom::area(X)
  
  if(is.null(closePairs)) {
    # Calculate pairwise relationships between all cells for a r
    closePairs <- spatstat.geom::closepairs(X, r, what = "ijd", distinct = FALSE) |> 
      data.frame()
  } else {
    closePairs <- closePairs |> 
      filter(d < r)
  }
  
  
  # Add cell type names to closePairs
  cellTypes <- data$cellType
  names(cellTypes) <- data$cellID
  closePairs$cellTypeI <- cellTypes[(closePairs$i)]
  closePairs$cellTypeJ <- cellTypes[(closePairs$j)]
  closePairs$i <- factor(closePairs$i, levels = data$cellID)

  # Count the number of each cell type near each cell.  
  # data.table would make this faster too.
  counts <- closePairs |>
    group_by(i, cellTypeI, cellTypeJ) |>
    summarise(n = sum(edge), .groups = "drop") |>
    pivot_wider(id_cols = c("i", "cellTypeI"), names_from = cellTypeJ, values_from = n,  values_fill = 0) |>
    as.data.frame() 
  
  ########
  # Calculate statistics
  ########
  
  Kontextual <- .Kontext(closePairs, counts, child1, child2, parent, r, Area, returnWeight)

  if (inhom) {
    L <- .Linhomfunction(closePairs, counts, child1, child2, r, Area) 
  } else{
    L <- .Lfunction(closePairs, counts, child1, child2, r, Area) 
  }
  
  if(returnWeight){
    return(Kontext)
  }
  
  return(c(L = L, Kontextual = Kontextual))
  
}



#' Calculates l function
#' @noRd
.Lfunction <- function(closePairs, counts, child1, child2, r, Area){
  
  nChild1 = sum(counts$cellTypeI == child1)
  nChild2 = sum(counts$cellTypeI == child2)
  
  # Adds up number of cell type 2 within cell type
  numerator <- sum(counts[counts$cellType==child1, child2])
  
  lambda2 = nChild2/Area
  kvalue = (numerator/lambda2) * (1/nChild1)
  
  sqrt(kvalue/pi) - r
  
}

#' Calculates inhomogenous L function
#' @noRd
.Linhomfunction <- function(closePairs, counts, child1, child2, r, Area){

  nChild1 = sum(counts$cellTypeI == child1)
  nChild2 = sum(counts$cellTypeI == child2)

  # Local intensity of each point
  weightChild1 = pull(counts, child1)/(pi*r^2)
  weightChild2 = pull(counts, child2)/(pi*r^2)

  closePairsShort = closePairs[closePairs$cellTypeI == child1 & closePairs$cellTypeJ == child2, ]
  closePairsShort$weightParent = closePairsShort$edge/(weightChild1[closePairsShort$i]* weightChild2[closePairsShort$j])

  numerator = sum(closePairsShort$weightParent)

  denominator = sum(1/weightChild1[counts$cellTypeI == child1])

  sqrt((numerator/denominator)/pi) - r

}




#' Calculates kontextual value
#' @noRd
.Kontext <- function(closePairs, counts, child1, child2, parent, r, Area, returnWeight){
  #child1 is root cell
  #child2 is child cell

  
  nParent = sum(counts$cellTypeI %in% parent)
  nChild1 = sum(counts$cellTypeI == child1)
  nChild2 = sum(counts$cellTypeI == child2)
  
  
  countsParent = counts |>
    select(any_of(parent)) |>
    rowSums()
  
  lambdaParent = countsParent/(pi*r^2)
  
  lambdaChild1 = lambdaParent
  lambdaChild2 = (lambdaParent/nParent)*nChild2
  
  closePairsShort <- closePairs[closePairs$cellTypeI==child1 & closePairs$cellTypeJ==child2,]

  closePairsShort$weightParent <- (closePairsShort$edge * lambdaChild1[closePairsShort$i])/(lambdaChild2[closePairsShort$j])
  
  
  numerator <- sum(closePairsShort$weightParent)
  
  # Weighted mean over child 2
  denominator <- sum(lambdaChild1[counts$cellTypeI == child1])
  
  # Turn this into centered L.
  sqrt(numerator/denominator/pi) - r
  
}





#' Converts a PPP object to a data frame for inhomLParent function
#' @noRd
#' @import spatstat.geom
#'
PPPdf <- function(ppp) {
  x <- as.data.frame(ppp)
  x$cellType <- factor(x$marks)
  x$x <- x$x
  x$y <- x$y
  x$cellID <- factor(seq_len(nrow(x)))
  return(x)
}


#' Edge correction for L values
#' @noRd
#' @import spatstat.geom

.borderEdge <- function(X, maxD) {
  W <- X$window
  bW <- spatstat.geom::union.owin(
    spatstat.geom::border(W, maxD, outside = FALSE),
    spatstat.geom::border(W, 2, outside = TRUE)
  )
  inB <- spatstat.geom::inside.owin(X$x, X$y, bW)
  e <- rep(1, X$n)
  if (any(inB)) {
    circs <- spatstat.geom::discs(X[inB], maxD, separate = TRUE)
    circs <-
      spatstat.geom::solapply(circs, spatstat.geom::intersect.owin, X$window)
    areas <- unlist(lapply(circs, spatstat.geom::area)) / (pi * maxD^2)
    e[inB] <- areas
  }

  e
}
