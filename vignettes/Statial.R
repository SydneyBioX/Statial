params <-
list(test = FALSE)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE
)
library(BiocStyle)

## ----warning = FALSE, message = FALSE-----------------------------------------
# Loading required packages
library(Statial)
library(spicyR)
library(ClassifyR)
library(lisaClust)
library(dplyr)
library(SingleCellExperiment)
library(ggplot2)
library(ggsurvfit)
library(survival)

theme_set(theme_classic())
nCores <- 4

## ----eval = FALSE-------------------------------------------------------------
#  # Install the package from Bioconductor
#  if (!requireNamespace("BiocManager", quietly = TRUE)) {
#    install.packages("BiocManager")
#  }
#  BiocManager::install("Statial")

