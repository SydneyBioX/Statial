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
library(spatstat)
library(EBImage)

theme_set(theme_classic())

## ----eval = FALSE-------------------------------------------------------------
#  # Install the development version from GitHub:
#  # install.packages("devtools")
#  devtools::install_github("SydneyBioX/Statial")

## -----------------------------------------------------------------------------
# Load head and neck data
sce = readRDS("../data/Head and Neck/sce.RDS")

# Examine all cell types in image
unique(sce$cellType)

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
konditionalResult = Konditional(
    imageData = sce,
    r = 50,
    from = "TC_CD4",
    to = "SC1",
    parent = immune,
    cores = 40
)

konditionalResult


## -----------------------------------------------------------------------------
# Get all relationships between cell types and their parents
parentDf = parentCombinations(all = all, tumour, bcells, tcells, myeloid, endothelial, epithelial, tissue, immune)
parentDf

## -----------------------------------------------------------------------------
# Selecting Image 1 as an example
image_1 = sce %>% colData %>% data.frame() %>% filter(imageID == "1")

konditionalResult = Konditional(image_1, 
                                parentDf = parentDf,
                                r = 50,
                                cores = 40)
konditionalResult

