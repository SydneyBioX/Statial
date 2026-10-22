# Plotting the original and kontextual L values over a range of radii.

This function takes outputs from rsCurve and plots them in ggplot. If
standard deviation is estimated in rsCurve, then confidence intervals
will be constructed based on the standard deviation. If the confidence
interval overlaps with 0, then the relationship is insignificant for
that radius.

## Usage

``` r
kontextPlot(rsDf)
```

## Arguments

- rsDf:

  A data frame from
  [`kontextCurve`](https://sydneybiox.github.io/Statial/reference/kontextCurve.md).

## Value

A ggplotly object showing the original and kontextual L function values
over a range of radii

## Examples

``` r
data("kerenSCE")

kerenImage6 <- kerenSCE[, kerenSCE$imageID == "6"]

rsDf <- kontextCurve(
  cells = kerenImage6,
  from = "p53",
  to = "Immune",
  parent = c("p53", "Keratin+Tumour"),
  rs = seq(10, 510, 100),
  cores = 2
)

kontextPlot(rsDf)
#> Ignoring unknown labels:
#> • fill : "Function"

```
