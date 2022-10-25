#' Creates a window for a PPP object
#'
#' @description
#' This function creates a window for a `spatstat::ppp` object, the type of 
#' window can be specified using the `window` argument.
#'
#' @param data A single image data frame from a SingleCellExperiment object or PPP object.
#' @param window The shape of window around the regions, can be `square`, `convex` or `concave`
#' @param window.length A tuning parameter for controlling the level of concavity when estimating concave windows.
#' @return Creates an `owin` class, representing the observation window for the image.
#'
#' @examples
#' data <- data.frame(x = rnorm(10), y = rnorm(10))
#' ow <- makeWindow(data, window = "square")
#'
#' spatstat.geom::ppp(x = data$x, y = data$y, window = ow)
#'
#' @export
#' @rdname makeWindow
#' @import spatstat
#' @import concaveman

makeWindow <- function(data,
                       window = "square",
                       window.length = NULL) {
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
    } else {
      window.length <- (max(data$x) - min(data$x)) / 20 * window.length
    }
    dist <- (max(data$x) - min(data$x)) / (length(data$x))
    bigDat <-
      do.call("rbind", lapply(as.list(as.data.frame(t(data[, c("x", "y")]))), function(x) {
        cbind(
          x[1] + c(0, 1, 0, -1, -1, 0, 1, -1, 1) * dist,
          x[2] + c(0, 1, 1, 1, -1, -1, -1, 0, 0) * dist
        )
      }))
    ch <-
      concaveman::concaveman(bigDat,
        length_threshold = window.length,
        concavity = 1
      )
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
