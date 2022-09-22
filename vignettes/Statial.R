params <-
list(test = FALSE)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(BiocStyle)

## ----warning=FALSE, message=FALSE---------------------------------------------
# load required packages
library(Statial)

# .libPaths("/dora/nobackup/biostat/Rpackages/v4")

library(ggplot2)
library(spicyR)
library(tidyverse)
library(SingleCellExperiment)

## ---- eval = FALSE------------------------------------------------------------
#  # Install the development version from GitHub:
#  # install.packages("devtools")
#  devtools::install_github("SydneyBioX/Statial")

## -----------------------------------------------------------------------------
# sce = readRDS("/albona/nobackup2/biostat/datasets/spatial/Damond2019_Diabetes_IMC/analysis/DamondCellsSCE.rds")

# cellExp = readRDS("../data/cellExp.RDS")
# intensities = data.frame(t(data.frame(cellMarks(cellExp))))
# colnames(intensities) =  paste0("cellID", cellSummary(cellExp)$cellID)
# 
# spatialMeasurements = data.frame(cellSummary(cellExp))
# rownames(spatialMeasurements) = paste0("cellID", cellSummary(cellExp)$cellID)
#   
# sce = SingleCellExperiment(list(intensities = intensities),
#                            colData = spatialMeasurements)
# 
# 
# saveRDS(sce, "../data/sce.RDS")


## -----------------------------------------------------------------------------
sce = readRDS("../data/sce.RDS")

## -----------------------------------------------------------------------------
intensitiesData = data.frame(t(assay(sce, "intensities")))
spatialData = data.frame(colData(sce))
markersToUse = colnames(intensitiesData)

singleCellData = cbind(spatialData[rownames(intensitiesData),], intensitiesData)
singleCellData = singleCellData %>% 
    mutate_at(markersToUse, function(x) ifelse(is.na(x), 0, x)) %>% 
    mutate_if(is.factor, as.character)

