#' Evaluation of Kontextual over a range of radii.
#' 
#' @description 
#' This function obtains `Kondtional` values over a range of radii, standard 
#' deviations for each value can be obtained using permutation for significance 
#' testing. To obtain estimates for standard deviations specify `se = TRUE`.
#'
#'
#' @param cells A single image from a SingleCellExperiment object
#' @param from The first cell type to be evaluated in the pairwise relationship.
#' @param to The second cell type to be evaluated in the pairwise relationship.
#' @param parent The parent population of the from cell type (must include from cell type).
#' @param rs A vector of radii to evaluate kontextual over.
#' @param inhom A logical value indicating whether to perform an inhomogeneous L function.
#' @param edge A logical value indicating whether to perform edge correction.
#' @param se A logical value to indicate if the standard deviation of
#' kontextual should be calculated to construct error bars.
#' @param nSim Number of randomisations to perform using \code{\link[Statial]{relabelKontextual}},
#'  which will be used to calculated the SE.
#' @param cores Number of cores for parallel processing.
#' @param ... Any arguments passed into \code{\link[Statial]{Kontextual}}.
#'
#' @return A data frame of original L values and Kontextual values evaluated
#' over a range of radii.
#'
#' @examples
#'
#' data("kerenSCE")
#' 
#' kerenImage6 = kerenSCE[, kerenSCE$imageID =="6"]
#' 
#' rsDf <- rsCurve(
#'   cells = kerenSCE,
#'   from = "CD4_Cell",
#'   to = "Keratin_Tumour",
#'   parent = c("CD4_Cell", "Macrophages"),
#'   rs = seq(10, 510, 100),
#'   cores = 2
#' )
#'
#' @export
#' @rdname rsCurve
#' @import tidyverse
#' @importFrom stats sd


rsCurve <- function(cells,
                    from,
                    to,
                    parent,
                    rs = seq(10, 100, 10),
                    inhom = TRUE,
                    edge = FALSE,
                    se = FALSE,
                    nSim = 20,
                    cores = 1,
                    ...) {
  kontextualVals <- Kontextual(
    cells = cells,
    from = from,
    to = to,
    parent = parent,
    r = rs,
    inhom = inhom,
    edgeCorrect = edge,
    cores = cores,
    includeOriginal = TRUE,
    ...
  )

  rsDf <- kontextualVals |>
    select("r", "original", "kontextual")


  if (se == TRUE) {
    seDf <- relabelKontextual(
      cells = cells,
      nSim = nSim,
      r = rs,
      from = from,
      to = to,
      parent = parent,
      returnImages = FALSE,
      inhom = inhom,
      edge = edge,
      cores = cores,
      ...
    )

    seDf <- seDf |>
      filter(type != "original") |>
      select("r", "original", "kontextual") |>
      group_by(r) |>
      summarise(
        "originalSd" = sd(original),
        "kontextualSd" = sd(kontextual)
      )

    rsDf <- merge(rsDf, seDf, by = "r")
  }

  return(rsDf)
}



#' Plotting the original and kontextual L values over a range of radii.
#'
#' @description 
#' This function takes outputs from \code{\link[Statial]{rsCurve}} and plots
#' them in ggplot. If standard deviation is estimated in \code{\link[Statial]{rsCurve}},
#' then confidence intervals will be constructed based on the standard deviation.
#' If the confidence interval overlaps with 0, then the relationship is insignificant 
#' for that radius.
#'
#' @param rsDf A data frame from \code{\link[Statial]{rsCurve}}.
#'
#' @return A ggplotly object showing the original and kontextual L function
#'  values over a range of radii
#'
#' @examples
#' data("kerenImage")
#'
#' rsDf <- rsCurve(
#'   cells = kerenImage,
#'   from = "p53",
#'   to = "Immune",
#'   parent = c("p53", "Keratin+Tumour"),
#'   rs = seq(10, 510, 100),
#'   cores = 2
#' )
#'
#' ggplotRs(rsDf)
#'
#' @export
#' @rdname ggplotRs
#' @import tidyverse
#' @import ggplot2
#' @importFrom stringr str_detect
#' @importFrom dplyr select
#' @importFrom tidyselect starts_with

ggplotRs <- function(rsDf) {
  if (str_detect(names(rsDf), "Sd") |> any()) {
    kontextual <- rsDf |>
      select(r, starts_with("kontextual")) |>
      mutate(
        "lower" = kontextual - kontextualSd,
        "upper" = kontextual + kontextualSd
      ) |>
      select(-"kontextualSd") |>
      pivot_longer(kontextual)


    original <- rsDf |>
      select("r", starts_with("original")) |>
      mutate(
        "lower" = original - originalSd,
        "upper" = original + originalSd
      ) |>
      select(-"originalSd") |>
      pivot_longer(original)

    seDf <- rbind(kontextual, original)


    p <- ggplot(seDf, aes(x = r, y = value, col = name)) +
      geom_point() +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = name, alpha = 0.2)) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "red"
      ) +
      geom_smooth(formula = y ~ x, method = "loess", se = FALSE) +
      guides(alpha = "none") +
      labs(
        x = "Radius (r)",
        y = "L(r) - r",
        fill = "Function",
        col = "Function"
      )
  } else {
    p <- rsDf |>
      pivot_longer(-r) |>
      ggplot(aes(x = r, y = value, col = name)) +
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_smooth(formula = y ~ x, method = "loess", se = FALSE) +
      labs(
        x = "Radius (r)",
        y = "L(r) - r",
        fill = "Function",
        col = "Function"
      )
  }

  return(p)
}
