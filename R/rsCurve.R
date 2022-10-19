#' Evaluation of Konditional over a range of radii.
#'
#'
#' @param cells A single image data frame from a SingleCellExperiment object
#' @param from The first cell type to be evaluated in the pairwise relationship.
#' @param to The second cell type to be evaluated in the pairwise relationship.
#' @param parent The parent population of the from cell type (must include from cell type).
#' @param rs A vector of radii to evaluate konditional over.
#' @param inhom A logical value indicating whether to perform an inhomogeneous L function.
#' @param edge A logical value indicating whether to perform edge correction.
#' @param se A logical value to indicate if the standard deviation of
#' konditional should be calculated to construct error bars.
#' @param nSim Number of randomisations to perform using \code{\link[Statial]{relabelKonditional}},
#'  which will be used to calculated the SE.
#' @param cores Number of cores for parallel processing.
#' @param ... Any arguments passed into \code{\link[Statial]{Konditional}}.
#'
#' @return A data frame of original L values and Konditional values evaluated
#' over a range of radii.
#'
#' @examples
#'
#' data("kerenImage")
#'
#' rsDf <- rsCurve(
#'   cells = kerenImage,
#'   from = "p53",
#'   to = "Immune",
#'   parent = c("p53", "Keratin+Tumour"),
#'   rs = seq(10, 510, 100),
#'   cores = 40
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
  konditionalVals <- Konditional(
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

  rsDf <- konditionalVals %>%
    select('r', 'original', 'konditional')


  if (se == TRUE) {
    seDf <- relabelKonditional(
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

    seDf <- seDf %>%
      filter(type != "original") %>%
      select('r', 'original', 'konditional') %>%
      group_by(r) %>%
      summarise(
        'originalSd' = sd(original),
        'konditionalSd' = sd(konditional)
      )

    rsDf <- merge(rsDf, seDf, by = "r")
  }

  return(rsDf)
}



#' Plotting the original and konditional L values over a range of radii.
#'
#'
#' @param rsDf A data frame from \code{\link[Statial]{rsCurve}}.
#'
#' @return A ggplotly object showing the original and konditional L function
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
#'   cores = 40
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
  if (str_detect(names(rsDf), "Sd") %>% any()) {
    konditional <- rsDf %>%
      select(r, starts_with("konditional")) %>%
      mutate(
        'lower' = konditional - konditionalSd,
        'upper' = konditional + konditionalSd
      ) %>%
      select(-'konditionalSd') %>%
      pivot_longer(konditional)


    original <- rsDf %>%
      select('r', starts_with("original")) %>%
      mutate(
        'lower' = original - originalSd,
        'upper' = original + originalSd
      ) %>%
      select(-'originalSd') %>%
      pivot_longer(original)

    seDf <- rbind(konditional, original)


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
        x = "Raidus (r)",
        y = "L(r) - r",
        fill = "Function",
        col = "Function"
      )
  } else {
    p <- rsDf %>%
      pivot_longer(-r) %>%
      ggplot(aes(x = r, y = value, col = name)) +
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_smooth(formula = y ~ x, method = "loess", se = FALSE) +
      labs(
        x = "Raidus (r)",
        y = "L(r) - r",
        fill = "Function",
        col = "Function"
      )
  }

  return(p)
}
