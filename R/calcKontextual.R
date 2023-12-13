#' Calculates l function
#' @noRd
.Lfunction <- function(closePairs, counts, child1, child2, r, area) {
  nChild1 <- sum(counts$cellTypeI == child1)
  nChild2 <- sum(counts$cellTypeI == child2)

  # Adds up number of cell type 2 within cell type
  numerator <- sum(counts[counts$cellType == child1, child2])

  lambda2 <- nChild2 / area
  kvalue <- (numerator / lambda2) * (1 / nChild1)

  sqrt(kvalue / pi) - r
}

#' Calculates inhomogenous L function
#' @noRd
.Linhomfunction <- function(closePairs, counts, child1, child2, r, area) {
  nChild1 <- sum(counts$cellTypeI == child1)
  nChild2 <- sum(counts$cellTypeI == child2)

  # Local intensity of each point
  weightChild1 <- pull(counts, child1) / (pi * r^2)
  weightChild2 <- pull(counts, child2) / (pi * r^2)

  closePairsShort <- closePairs[closePairs$cellTypeI == child1 & closePairs$cellTypeJ == child2, ]
  closePairsShort$weightParent <- closePairsShort$edge / (weightChild1[closePairsShort$i] * weightChild2[closePairsShort$j])

  numerator <- sum(closePairsShort$weightParent)

  denominator <- sum(1 / weightChild1[counts$cellTypeI == child1])

  sqrt((numerator / denominator) / pi) - r
}




#' Calculates kontextual value
#' @noRd
.Kontext <- function(closePairs, counts, child1, child2, parent, r, returnWeight) {
  # child1 is root cell
  # child2 is child cell


  nParent <- sum(counts$cellTypeI %in% parent)
  nChild1 <- sum(counts$cellTypeI == child1)
  nChild2 <- sum(counts$cellTypeI == child2)


  countsParent <- counts |>
    select(any_of(parent)) |>
    rowSums()

  lambdaParent <- countsParent / (pi * r^2)

  lambdaChild1 <- lambdaParent
  lambdaChild2 <- (lambdaParent / nParent) * nChild2

  closePairsShort <- closePairs[closePairs$cellTypeI == child1 & closePairs$cellTypeJ == child2, ]

  closePairsShort$weightParent <- (closePairsShort$edge * lambdaChild1[closePairsShort$i]) / (lambdaChild2[closePairsShort$j])


  numerator <- sum(closePairsShort$weightParent)

  # Weighted mean over child 2
  denominator <- sum(lambdaChild1[counts$cellTypeI == child1])

  # Turn this into centered L.
  sqrt(numerator / denominator / pi) - r
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
