# Extract parent and all children from a Phylo object

This function takes in a \`phylo\` object or a \`treekoR\` result from
the
[`getClusterTree`](https://rdrr.io/pkg/treekoR/man/getClusterTree.html)
function, and converts its into a named list of each and children to
input into
[`parentCombinations`](https://sydneybiox.github.io/Statial/reference/parentCombinations.md).

**Note**: Parent populations with one child will be pruned. Make sure to
include this cell type in the \`all\` vector when using
[`parentCombinations`](https://sydneybiox.github.io/Statial/reference/parentCombinations.md)
to ensure this cell type is included in pairwise calculations.

## Usage

``` r
getParentPhylo(phylo_tree)
```

## Arguments

- phylo_tree:

  a phylo object or a treekoR result.

## Value

A named list of parents and their respective children.
