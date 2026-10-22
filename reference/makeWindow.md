# Creates a window for a PPP object

This function creates a window for a \`spatstat::ppp\` object, the type
of window can be specified using the \`window\` argument.

## Usage

``` r
makeWindow(data, window = "square", window.length = NULL)
```

## Arguments

- data:

  A single image data frame from a SingleCellExperiment object or PPP
  object.

- window:

  The shape of window around the regions, can be \`square\`, \`convex\`
  or \`concave\`

- window.length:

  A tuning parameter for controlling the level of concavity when
  estimating concave windows.

## Value

Creates an \`owin\` class, representing the observation window for the
image.

## Examples

``` r
data <- data.frame(x = rnorm(10), y = rnorm(10))
ow <- makeWindow(data, window = "square")

spatstat.geom::ppp(x = data$x, y = data$y, window = ow)
#> Planar point pattern: 10 points
#> window: rectangle = [-1.8218177, 2.0650249] x [-1.8630115, 1.8885049] units
```
