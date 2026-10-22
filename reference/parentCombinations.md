# Create all combinations of cell type relationships from a list of parents

This function takes in named vectors of all the parent populations in
the dataset, and creates a data frame containing all pairwise cell
relationships, this data frame can be inputed into the \`parentDf\`
argument in \`Kontextual\`.

## Usage

``` r
parentCombinations(all, ..., parentList = NULL)
```

## Arguments

- all:

  A list of all the \`to\` cell types Kontextual is evaluated over.

- ...:

  Vectors of each parent population.

- parentList:

  a named list where the names correspond to parent names and the values
  contain a vector of children for that parent. Note: If parentList is
  specified the \`...\` argument will be ignored, see examples.

## Value

A data frame containing all pairwise cell relationships and their
corresponding parent

## Examples

``` r
# Example 1, using `parentList`

parentList <- list(
  "tcells" = c("CD4", "CD8"),
  "tissue" = c("epithelial", "stromal")
)

allCells <- c("tumour", "CD4", "CD8", "epithelial", "stromal")

parentCombinations(all = allCells, parentList = parentList)
#>          from         to       parent parent_name
#> 1      tumour        CD4     CD4, CD8      tcells
#> 2         CD8        CD4     CD4, CD8      tcells
#> 3  epithelial        CD4     CD4, CD8      tcells
#> 4     stromal        CD4     CD4, CD8      tcells
#> 5      tumour        CD8     CD4, CD8      tcells
#> 6         CD4        CD8     CD4, CD8      tcells
#> 7  epithelial        CD8     CD4, CD8      tcells
#> 8     stromal        CD8     CD4, CD8      tcells
#> 9      tumour epithelial epitheli....      tissue
#> 10        CD4 epithelial epitheli....      tissue
#> 11        CD8 epithelial epitheli....      tissue
#> 12    stromal epithelial epitheli....      tissue
#> 13     tumour    stromal epitheli....      tissue
#> 14        CD4    stromal epitheli....      tissue
#> 15        CD8    stromal epitheli....      tissue
#> 16 epithelial    stromal epitheli....      tissue


# Example 2, with `...` operator
tcells <- c("CD4", "CD8")
tissue <- c("epithelial", "stromal")
allCells <- c("tumour", tissue, tcells)

parentCombinations(all = allCells, tcells, tissue)
#>          from         to       parent parent_name
#> 1      tumour        CD4     CD4, CD8      tcells
#> 2  epithelial        CD4     CD4, CD8      tcells
#> 3     stromal        CD4     CD4, CD8      tcells
#> 4         CD8        CD4     CD4, CD8      tcells
#> 5      tumour        CD8     CD4, CD8      tcells
#> 6  epithelial        CD8     CD4, CD8      tcells
#> 7     stromal        CD8     CD4, CD8      tcells
#> 8         CD4        CD8     CD4, CD8      tcells
#> 9      tumour epithelial epitheli....      tissue
#> 10    stromal epithelial epitheli....      tissue
#> 11        CD4 epithelial epitheli....      tissue
#> 12        CD8 epithelial epitheli....      tissue
#> 13     tumour    stromal epitheli....      tissue
#> 14 epithelial    stromal epitheli....      tissue
#> 15        CD4    stromal epitheli....      tissue
#> 16        CD8    stromal epitheli....      tissue
```
