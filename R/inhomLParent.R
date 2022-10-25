#' @noRd
#' @import spatstat
#' @importFrom methods is
#' @importFrom stats quantile
inhomLParent <- function(data,
                         Rs = 20,
                         window = "convex",
                         window.length = NULL,
                         weightQuantile = .80,
                         from = NULL,
                         to = NULL,
                         parent = NULL,
                         inhom = TRUE,
                         edgeCorrect = TRUE,
                         includeZeroCells = TRUE,
                         closePairs = NULL) {
  if (is(data, "ppp")) {
    data <- PPPdf(data)
  }

  if (!("cellID" %in% names(data))) {
    data$cellID <- factor(seq_len(nrow(data)))
  }
  # if class is data frame it show make window etc.
  ow <- makeWindow(data, window, window.length)

  # sigma = radius
  if (is.null(Rs)) {
    Rs <- 50
  }
  sigma <- Rs

  X <-
    spatstat.geom::ppp(
      x = data$x,
      y = data$y,
      window = ow,
      marks = data$cellType
    )


  Y <- X
  if (!is.null(parent)) Y <- Y[Y$marks %in% parent, ]

  Area <- area(X)


  if (!is.null(Rs)) {
    # Use disc kernel for interpretation
    den <-
      spatstat.core::density.ppp(Y, sigma = sigma, kernel = "disc")
    den <- den / max(den)
    # den$v <- pmax(den$v, minLambda)
    Area <- area(X) * mean(den) # area of the parent
  }


  maxR <- min(ow$xrange[2] - ow$xrange[1], ow$yrange[2] - ow$yrange[1]) / 2.01
  Rs <- unique(pmin(c(0, sort(Rs)), maxR))


  if (is.null(from)) from <- levels(data$cellType)
  if (is.null(to)) to <- levels(data$cellType)

  use <- data$cellType %in% c(from, to)
  fulldata <- data
  data <- data[use, ]
  X <- X[use, ]

  # pairwise relationships for a r
  if (is.null(closePairs)) {
    closePairs <-
      spatstat.geom::closepairs(X, max(Rs), what = "ijd", distinct = FALSE)
    closePairs$j <- data$cellID[closePairs$j]
    closePairs$i <- data$cellID[closePairs$i]
  } else {
    closePairs <-
      closePairs |>
      data.frame() |>
      filter(d <= max(Rs))

    # Convert True false into indexes using seq_along and then subset closePairs dataframe

    closePairs <- closePairs[closePairs$i %in% seq_along(use)[use], ]
    closePairs <- closePairs[closePairs$j %in% seq_along(use)[use], ]
  }



  p <- closePairs

  n <- X$n

  cT <- data$cellType
  names(cT) <- data$cellID

  p$d <- cut(p$d, Rs, labels = Rs[-1], include.lowest = TRUE)

  w <- rep(1, length(X))
  p$wt <- rep(1, length(p$d))
  if (!is.null(sigma)) {
    np <- spatstat.geom::nearest.valid.pixel(X$x, X$y, den)
    w <- den$v[cbind(np$row, np$col)]
    names(w) <- data$cellID
    invWeight <- 1 / w[as.character(p$i)]
    p$wt <- pmin(invWeight, quantile(invWeight, weightQuantile))
  }

  num <- tapply(w, data$cellType, length)
  # weights of each cell type
  lam <- tapply(w, data$cellType, sum) / Area
  # sum of the weights
  if (inhom) {
    num <- tapply(pmin(1 / w, quantile(1 / w, weightQuantile)), data$cellType, sum)
  }
  if (!inhom) p$wt <- 1

  p$cellTypeJ <- cT[as.character(p$j)]
  p$cellTypeI <- cT[as.character(p$i)]
  p$i <- factor(p$i, levels = data$cellID)

  if (edgeCorrect) {
    edge <- borderEdge(X, Rs[-1])
    edge <- do.call("cbind", edge)
    edge <- as.data.frame(edge)
    colnames(edge) <- Rs[-1]
    edge$i <- data$cellID
    # edge$i <- factor(data$cellID, levels = data$cellID)
    edge <- tidyr::pivot_longer(edge, -i, "d")
    p <- dplyr::left_join(as.data.frame(p), edge, c("i", "d"))
  } else {
    p <- as.data.frame(p)
    p$value <- 1
  }


  p$d <- factor(p$d, levels = Rs[-1])

  p <- p[as.character(p$i) != as.character(p$j), ]

  use <- p$cellTypeI %in% from & p$cellTypeJ %in% to
  p <- p[use, ]

  r <- inhomL(p, lam, X, Rs, num, Area)

  wt <- r$wt
  names(wt) <- paste(r$cellTypeI, r$cellTypeJ, sep = "__")

  m1 <- rep(from, times = length(to))
  m2 <- rep(to, each = length(from))
  labels <- paste(m1, m2, sep = "__")

  assoc <- rep(-sum(Rs), length(labels))
  names(assoc) <- labels
  if (!includeZeroCells) assoc[!(m1 %in% X$marks & m2 %in% X$marks)] <- NA
  assoc[names(wt)] <- wt
  names(assoc) <- labels

  assoc
}



#' Calculates L value from weights and lambda values
#' @noRd
#' @importFrom data.table as.data.table setkey CJ .SD ":="

inhomL <-
  function(p, lam, X, Rs, num, Area) {
    r <- data.table::as.data.table(p)

    r$wt <-
      (r$wt) / as.numeric(lam[r$cellTypeJ]) / (as.numeric(num[r$cellTypeI]))

    r <- r[, j := NULL]
    r <- r[, value := NULL]
    r <- r[, i := NULL]
    data.table::setkey(r, d, cellTypeI, cellTypeJ)
    r <- r[data.table::CJ(d, cellTypeI, cellTypeJ, unique = TRUE)][, lapply(.SD, sum), by = .(d, cellTypeI, cellTypeJ)][is.na(wt), wt := 0]
    r <- r[, wt := cumsum(wt), by = list(cellTypeI, cellTypeJ)]
    r <- r[, list(wt = sum(sqrt(wt / pi))), by = .(cellTypeI, cellTypeJ)]
    r$wt <- r$wt - sum(Rs)

    r <- as.data.frame(r)

    r
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

borderEdge <- function(X, maxD) {
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
