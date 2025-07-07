#' Create all combinations of cell type relationships from a list of parents
#'
#' @description
#' This function takes in named vectors of all the parent populations in the
#' dataset, and creates a data frame containing all pairwise cell relationships,
#' this data frame can be inputed into the `parentDf` argument in `Kontextual`.
#'
#' @param all A list of all the `to` cell types Kontextual is evaluated over.
#' @param parentList a named list where the names correspond to parent names and
#' the values contain a vector of children for that parent. Note: If parentList
#' is specified the `...` argument will be ignored, see examples.
#' @param ... Vectors of each parent population.
#'
#' @return A data frame containing all pairwise cell relationships and their
#' corresponding parent
#'
#' @examples
#' # Example 1, using `parentList`
#'
#' parentList <- list(
#'   "tcells" = c("CD4", "CD8"),
#'   "tissue" = c("epithelial", "stromal")
#' )
#' 
#' allCells <- c("tumour", "CD4", "CD8", "epithelial", "stromal")
#' 
#' parentCombinations(all = allCells, parentList = parentList)
#' 
#' 
#' # Example 2, with `...` operator
#' tcells <- c("CD4", "CD8")
#' tissue <- c("epithelial", "stromal")
#' allCells <- c("tumour", tissue, tcells)
#'
#' parentCombinations(all = allCells, tcells, tissue)
#'
#' @export
#' @rdname parentCombinations
#' @import dplyr
#' @import tidyr
#'
parentCombinations <- function(all, ..., parentList = NULL) {
  
  if(is.null(parentList)) {
    # Gets variable names of all the parent vector
    names <- as.list(substitute(c(...)))[-1]
    
    parentList <- list(...)
    names(parentList) <- names
  }

  # Creates data.frame of parent name and parent vector
  parentTable <- data.frame(
    parent_name = names(parentList),
    parent = I(parentList)
  )


  # Creates all combination of parent and child
  parentDfs <- lapply(
    seq_along(parentList),
    function(x) {
      return(tidyr::crossing(
        to = parentList[[x]],
        parent_name = names(parentList)[x]
      ))
    }
  )

  parentDf <- bind_rows(parentDfs) |>
    merge(parentTable, by = "parent_name") |>
    expand_grid(from = unique(all)) |>
    data.frame() |>
    select("from", "to", "parent", "parent_name") |>
    filter(from != to)



  return(parentDf)
}




#' Extract parent and all children from a Phylo object
#'
#' @description
#' This function takes in a `phylo` object or a `treekoR` result from the 
#' \code{\link[treekoR]{getClusterTree}} function, and converts its into a 
#' named list of each and children to input into 
#' \code{\link[Statial]{parentCombinations}}. 
#' 
#' \bold{Note}: Parent populations with one child will be pruned. Make sure to 
#' include this cell type in the `all` vector when using 
#' \code{\link[Statial]{parentCombinations}} to ensure this cell type is included
#' in pairwise calculations.
#
#'
#' @param phylo_tree a phylo object or a treekoR result.
#'
#' @return A named list of parents and their respective children.
#'
#'
#' @export
#' @rdname getParentPhylo
#' @importFrom dplyr rename filter mutate group_by summarise
getParentPhylo = function(phylo_tree) {
  
  # Get tree from treekor object
  if(class(phylo_tree) == "list" && "clust_tree" %in% names(phylo_tree)) {
    phylo_tree = phylo_tree$clust_tree
  }
  
  if(!class(phylo_tree) == "phylo") {
    stop("Please input a phylo object or a treekoR result")
  }
  
  node_labels <- c(phylo_tree$tip.label)
  node_length = length(node_labels)
  names(node_labels) = 1:node_length
  
  edge_matrix <- data.frame(phylo_tree$edge) |> 
    rename("parent" = "X1",
           "child" = "X2") |> 
    filter(child <= node_length) |> 
    mutate(child = node_labels[child]) |>
    mutate(child = unname(child)) |> 
    group_by(parent) |> 
    summarise(children = list(child)) |> 
    filter(lengths(children) > 1)
  
  child_list = edge_matrix$children
  
  names(child_list) = paste("parent", seq_along(child_list), sep = "_")
  
  return(child_list)
}