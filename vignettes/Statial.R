params <-
list(test = FALSE)

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE,
  message = FALSE,
  warning = FALSE
)
library(BiocStyle)

## ----warning = FALSE, message = FALSE, eval = TRUE----------------------------
# Loading required packages
library(Statial)
#library(spicyR)
devtools::load_all("/dski/nobackup/farhana/spicyR")
library(tidyverse)
library(ggthemes)
library(ggh4x)

library(ClassifyR)
library(lisaClust)
library(dplyr)
library(SingleCellExperiment)
library(ggplot2)
library(ggsurvfit)
library(survival)
library(tibble)
library(treekoR)

theme_set(theme_classic())
nCores <- 1

## -----------------------------------------------------------------------------
signifPlot(survivalResults)

## -----------------------------------------------------------------------------
# Extracting survival data
survData <- kerenSCE |>
  colData() |>
  data.frame() |>
  select(imageID, Survival_days_capped, Censored) |>
  mutate(event = 1 - Censored) |>
  unique()

# Creating survival vector
kerenSurv <- Surv(survData$Survival_days_capped, survData$event)
names(kerenSurv) <- survData$imageID

kerenSurv

## ----fig.width=5, fig.height=4------------------------------------------------
# Selecting most significant relationship
survRelationship <- kontextMat[["DC__NK__immune"]]
survRelationship <- ifelse(survRelationship > median(survRelationship), "Localised", "Dispersed")

# Plotting Kaplan-Meier curve
survfit2(kerenSurv ~ survRelationship) |>
  ggsurvfit() +
  ggtitle("DC__NK__immune")

## -----------------------------------------------------------------------------
kerenSCE <- getDistances(kerenSCE,
  maxDist = 200,
  nCores = 1
)

kerenSCE <- getAbundances(kerenSCE,
  r = 200,
  nCores = 1
)

## -----------------------------------------------------------------------------
stateChanges <- calcStateChanges(
  cells = kerenSCE,
  type = "distances",
  image = "6",
  from = "Keratin_Tumour",
  to = "Macrophages",
  marker = "p53",
  nCores = 1
)

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
  method = "lm"
)

p

## -----------------------------------------------------------------------------
stateChanges <- calcStateChanges(
  cells = kerenSCE,
  type = "distances",
  nCores = 1,
  minCells = 100
)

stateChanges |>
  filter(imageID == 6) |>
  head(n = 10)

## -----------------------------------------------------------------------------
p <- plotStateChanges(
  cells = kerenSCE,
  type = "distances",
  image = "6",
  from = "Keratin_Tumour",
  to = "Macrophages",
  marker = "HLA_Class_1",
  size = 1,
  shape = 19,
  interactive = FALSE,
  plotModelFit = FALSE,
  method = "lm"
)

p

## -----------------------------------------------------------------------------
stateChanges |> head(n = 10)

## -----------------------------------------------------------------------------
p <- plotStateChanges(
  cells = kerenSCE,
  type = "distances",
  image = "35",
  from = "CD4_Cell",
  to = "B",
  marker = "CD20",
  size = 1,
  shape = 19,
  interactive = FALSE,
  plotModelFit = FALSE,
  method = "lm"
)

p

## -----------------------------------------------------------------------------
kerenSCE <- calcContamination(kerenSCE)

stateChangesCorrected <- calcStateChanges(
  cells = kerenSCE,
  type = "distances",
  nCores = 1,
  minCells = 100,
  contamination = TRUE
)

stateChangesCorrected |> head(n = 20)

## -----------------------------------------------------------------------------
cellTypeMarkers <- c("CD3", "CD4", "CD8", "CD56", "CD11c", "CD68", "CD45", "CD20")

values <- c("blue", "red")
names(values) <- c("None", "Corrected")

df <- rbind(
  data.frame(TP = cumsum(stateChanges$marker %in% cellTypeMarkers), FP = cumsum(!stateChanges$marker %in% cellTypeMarkers), type = "None"),
  data.frame(TP = cumsum(stateChangesCorrected$marker %in% cellTypeMarkers), FP = cumsum(!stateChangesCorrected$marker %in% cellTypeMarkers), type = "Corrected")
)

ggplot(df, aes(x = TP, y = FP, colour = type)) +
  geom_line() +
  labs(y = "Cell state marker", x = "Cell type marker") +
  scale_colour_manual(values = values)

## -----------------------------------------------------------------------------
ggplot(df, aes(x = TP, y = FP, colour = type)) +
  geom_line() +
  xlim(0, 100) +
  ylim(0, 1000) +
  labs(y = "Cell state marker", x = "Cell type marker") +
  scale_colour_manual(values = values)

## -----------------------------------------------------------------------------
# Preparing features for Statial
stateMat <- prepMatrix(stateChanges)

# Ensuring rownames of stateMat match up with rownames of the survival vector
stateMat <- stateMat[names(kerenSurv), ]

# Remove some very small values
stateMat <- stateMat[, colMeans(abs(stateMat) > 0.0001) > .8]

survivalResults <- colTest(stateMat, kerenSurv, type = "survival")

head(survivalResults)

## ----fig.width=5, fig.height=4------------------------------------------------
# Selecting the most significant relationship
survRelationship <- stateMat[["Keratin_Tumour__CD4_Cell__Keratin6"]]
survRelationship <- ifelse(survRelationship > median(survRelationship), "Higher expression in close cells", "Lower expression in close cells")

# Plotting Kaplan-Meier curve
survfit2(kerenSurv ~ survRelationship) |>
  ggsurvfit() +
  add_pvalue() +
  ggtitle("Keratin_Tumour__CD4_Cell__Keratin6")

## -----------------------------------------------------------------------------
set.seed(51773)

# Preparing features for lisaClust
kerenSCE <- lisaClust::lisaClust(kerenSCE, k = 5)

## ----fig.height=5, fig.width=6.5----------------------------------------------
# Use hatching to visualise regions and cell types.
lisaClust::hatchingPlot(kerenSCE,
  useImages = "5",
  line.spacing = 41, # spacing of lines
  nbp = 100 # smoothness of lines
)

## ----lisaClust----------------------------------------------------------------
cellTypeRegionMeans <- getMarkerMeans(kerenSCE,
  imageID = "imageID",
  cellType = "cellType",
  region = "region"
)

survivalResults <- colTest(cellTypeRegionMeans[names(kerenSurv), ], kerenSurv, type = "survival")

head(survivalResults)

## -----------------------------------------------------------------------------
# Calculate cell type and region proportions
cellTypeProp <- getProp(kerenSCE,
  feature = "cellType",
  imageID = "imageID"
)
regionProp <- getProp(kerenSCE,
  feature = "region",
  imageID = "imageID"
)

# Combine all the features into a list for classification
featureList <- list(
  states = stateMat,
  kontextual = kontextMat,
  regionMarkerMeans = cellTypeRegionMeans,
  cellTypeProp = cellTypeProp,
  regionProp = regionProp
)

# Ensure the rownames of the features match the order of the survival vector
featureList <- lapply(featureList, function(x) x[names(kerenSurv), ])


set.seed(51773)

kerenCV <- crossValidate(
  measurements = featureList,
  outcome = kerenSurv,
  classifier = "CoxPH",
  selectionMethod = "CoxPH",
  nFolds = 5,
  nFeatures = 10,
  nRepeats = 20,
  nCores = 1
)

## -----------------------------------------------------------------------------
# Calculate AUC for each cross-validation repeat and plot.
performancePlot(kerenCV,
  characteristicsList = list(x = "Assay Name")
) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

## -----------------------------------------------------------------------------
sessionInfo()

