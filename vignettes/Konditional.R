params <-
list(test = FALSE)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(BiocStyle)

## ----warning = FALSE, message = FALSE-----------------------------------------
#Loading required packages
library(Statial)
library(tidyverse)
library(SingleCellExperiment)

theme_set(theme_classic())

## ----eval = FALSE-------------------------------------------------------------
#  # Install the development version from GitHub:
#  # install.packages("devtools")
#  devtools::install_github("SydneyBioX/Statial")
#  

## -----------------------------------------------------------------------------
# Load head and neck data
data("headSCE")

# Examine all cell types in image
unique(headSCE$cellType)

# Set up cell populations
tumour = c("SC1", "SC2", "SC3", "SC4", "SC5", "SC6", "SC7")

bcells = c("BC1", "BC2", "BC3")
tcells = c("TC_CD4", "TC_CD8")
myeloid = c("MC1", "MC2", "MC3")

endothelial = c("EC1", "EC2")
epithelial = c("EP")

tissue = c(endothelial, epithelial)
immune = c(bcells, tcells, myeloid, "GC") #GC = granulocytes

all = c(tumour, tissue, immune, "Undefined")



## -----------------------------------------------------------------------------
CD4_Konditional = Konditional(
    cells = headSCE,
    r = 50,
    from = "TC_CD4",
    to = "SC5",
    parent = immune,
    cores = 40
)

head(CD4_Konditional)

## ----fig.wide = TRUE----------------------------------------------------------

ggplot(CD4_Konditional, aes(x = original, y = konditional, col = imageID)) + 
    geom_point() +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    geom_vline(xintercept = 0, col = "red", linetype = "dashed")

## -----------------------------------------------------------------------------
# Get all relationships between cell types and their parents
parentDf = parentCombinations(all = all,
                              tumour,
                              bcells,
                              tcells,
                              myeloid,
                              endothelial,
                              epithelial,
                              tissue,
                              immune)
head(parentDf)

## -----------------------------------------------------------------------------
# Selecting Image 1 as an example
image_1 = headSCE %>% colData %>% data.frame() %>% filter(imageID == "1")


image1_Konditional = Konditional(cells = image_1,
                                 parentDf = parentDf,
                                 r = 50,
                                 cores = 40)
head(image1_Konditional)

## -----------------------------------------------------------------------------
data("kerenImage")

kerenImage %>% filter(cellType %in% c("Keratin+Tumour", "Immune", "p53")) %>%
    arrange(cellType) %>% 
  ggplot(aes(x = x, y = y, color = cellType)) + 
  geom_point(size = 1) +
  scale_colour_manual(values = c("#505050", "#D6D6D6", "#64BC46"))

## -----------------------------------------------------------------------------
rsDf = rsCurve(
    cells = kerenImage,
    from = "p53",
    to = "Immune",
    parent = c("p53", "Keratin+Tumour"),
    rs = seq(10, 510, 100),
    cores = 40
)


ggplotRs(rsDf)


## -----------------------------------------------------------------------------
sessionInfo()

