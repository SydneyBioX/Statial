# Convert Kontextual or state changes result to a matrix for classification

Convert Kontextual or state changes result to a matrix for
classification

## Usage

``` r
prepMatrix(result, replaceVal = 0, column = NULL, test = NULL)
```

## Arguments

- result:

  a kontextual or state changes result data.frame.

- replaceVal:

  value which NAs are replaced with.

- column:

  The column which contains the scores that you want to select.

- test:

  A column containing which will be the column names of the expanded
  matrix.

## Examples

``` r
data("kerenSCE")


CD4_Kontextual <- Kontextual(
  cells = kerenSCE,
  r = 50,
  from = "Macrophages",
  to = "Keratin_Tumour",
  parent = c("Macrophages", "CD4_Cell"),
  image = "6"
)


kontextMat <- prepMatrix(CD4_Kontextual)
```
