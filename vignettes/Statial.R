params <-
list(test = FALSE)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
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

## -----------------------------------------------------------------------------
# Load head and neck data
data("kerenSCE")

kerenSCE

## -----------------------------------------------------------------------------

# Examine all cell types in image
unique(kerenSCE$cellType)

# Set up cell populations
tumour <- c("Keratin_Tumour", "Tumour")

bcells <- c("B")
tcells <- c("CD3_Cell", "CD4_Cell", "CD8_Cell", "Tregs")
myeloid <- c("Dc/Mono", "DC", "Mono/Neu", "Macrophages", "other immune", "Neutrophils")

endothelial <- c("Endothelial")
mesenchymal <- c("Mesenchymal")

tissue <- c(endothelial, mesenchymal)
immune <- c(bcells, tcells, myeloid, "NK") # GC = granulocytes

all <- c(tumour, tissue, immune, "Unidentified")

## -----------------------------------------------------------------------------
#Select image 6 from the kerenSCE dataset
kerenImage6 = kerenSCE[, kerenSCE$imageID =="6"]

#Select for all cells that express higher than baseline level of p53
p53Pos = assay(kerenImage6)["p53",]  |> 
  as.numeric() > -0.300460
kerenImage6$cellType[p53Pos & kerenImage6$cellType %in% c("Keratin_Tumour")] <- "p53+Tumour"

#Group all immune cells under the name "Immune"
kerenImage6$cellType[kerenImage6$cellType %in% immune] <- "Immune"

kerenImage6 |>
  colData() %>%
  as.data.frame() %>%
  filter(cellType %in% c("Keratin_Tumour", "Immune", "p53+Tumour")) %>%
  arrange(cellType) %>%
  ggplot(aes(x = x, y = y, color = cellType)) +
  geom_point(size = 1) +
  scale_colour_manual(values = c("#505050", "#D6D6D6", "#64BC46"))


## -----------------------------------------------------------------------------

#Select for all cells that express higher than baseline level of p53
kerenSCE$cellTypeNew <- kerenSCE$cellType
p53Pos = assay(kerenSCE)["p53",]  |> 
  as.numeric() > -0.300460
kerenSCE$cellTypeNew[p53Pos & kerenSCE$cellType %in% c("Keratin_Tumour")] <- "p53+Tumour"

#Group all immune cells under the name "Immune"
kerenSCE$cellTypeNew[kerenSCE$cellType %in% immune] <- "Immune"


curves <- kontextCurve(
  cells = kerenSCE,
  from = "p53+Tumour",
  to = "Immune",
  parent = c("p53+Tumour", "Keratin_Tumour"),
  rs = seq(10, 510, 100),
  image = "6",
  cellType = "cellTypeNew",
  cores = nCores
)

kontextPlot(curves)

## -----------------------------------------------------------------------------
p53_Kontextual <- Kontextual(
  cells = kerenSCE,
  r = 50,
  from = "p53+Tumour",
  to = "Immune",
  parent = c("p53", "Keratin_Tumour"),
)

p53_Kontextual


## -----------------------------------------------------------------------------
# Get all relationships between cell types and their parents
parentDf <- parentCombinations(
  all = all,
  tumour,
  bcells,
  tcells,
  myeloid,
  endothelial,
  mesenchymal,
  tissue,
  immune
)


## ----eval=FALSE---------------------------------------------------------------
#  # Running Kontextual on all relationships across all images.
#  kerenKontextual <- Kontextual(
#    cells = kerenSCE,
#    parentDf = parentDf,
#    r = 50,
#    cores = nCores
#  )

## -----------------------------------------------------------------------------
data("kerenKontextual")
head(kerenKontextual)

## -----------------------------------------------------------------------------

kerenSCE <- getDistances(kerenSCE,
                    maxDist = 200,
                    nCores = 1)

kerenSCE <- getAbundances(kerenSCE,
                     r = 200,
                     nCores = 1)


## -----------------------------------------------------------------------------
#This function takes approximately 5 minutes to run with 40 cores. 
stateChanges <- calcStateChanges(
  cells = kerenSCE,
  type = "distances",
  image = "6",
  from = "Keratin_Tumour",
  to = "Macrophages",
  marker = "p53",
  nCores = 1)

stateChanges

## -----------------------------------------------------------------------------
p <- plotStateChanges(
  cells = kerenSCE,
  type = "distances",
  image = "6",
  from = "Keratin_Tumour",
  to = "Macrophages",
  marker = "p53",
  size = 1,
  shape = 19,
  interactive = FALSE,
  plotModelFit = FALSE,
  method = "lm")

p


## -----------------------------------------------------------------------------
#This function takes approximately 5 minutes to run with 40 cores. 
stateChanges <- calcStateChanges(
  cells = kerenSCE,
  type = "distances",
  nCores = 1,
  minCells = 100)

stateChanges |> head(n = 10)

stateChanges |> 
  filter(imageID == 6) |>
  head(n = 10)

## -----------------------------------------------------------------------------

p <- plotStateChanges(
  cells = kerenSCE,
  type = "distances",
  image = "27",
  from = "Keratin_Tumour",
  to = "Macrophages",
  marker = "HLA_Class_1",
  size = 1,
  shape = 19,
  interactive = FALSE,
  plotModelFit = FALSE,
  method = "lm")

p


