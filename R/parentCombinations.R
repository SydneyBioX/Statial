#' Create all combinations of cell type relationships from a list of parents
#'
#' @description
#' This function takes in named vectors of all the parent populations in the
#' dataset, and creates a data frame containing all pairwise cell relationships,
#' this data frame can be inputed into the `parentDf` argument in `Kontextual`.
#'
#' @param all A list of all the `to` cell types Kontextual is evaluated over
#' @param ... Vectors of each parent population
#'
#' @return A data frame containing all pairwise cell relationships and their
#' corresponding parent
#'
#' @examples
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
parentCombinations <- function(all, ...) {
  # Gets variable names of all the parent vector
  names <- as.list(substitute(c(...)))[-1]

  parentList <- list(...)
  names(parentList) <- names

  # Creates data.frame of parent name and parent vector
  parentTable <- data.frame(
    parent_name = names(parentList),
    parent = I(parentList)
  )


  # Creates all combination of parent and child
  parentDfs <- lapply(
    seq_along(parentList),
    function(x) {
      return(crossing(
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
