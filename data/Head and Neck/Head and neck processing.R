library(SingleCellExperiment)
cellExp = readRDS("data/Head and Neck/cellExp.RDS")
intensities = data.frame(t(data.frame(cellMarks(cellExp))))
colnames(intensities) =  paste0("cellID", cellSummary(cellExp)$cellID)

spatialMeasurements = data.frame(cellSummary(cellExp))
rownames(spatialMeasurements) = paste0("cellID", cellSummary(cellExp)$cellID)

sce = SingleCellExperiment(list(intensities = intensities),
                           colData = spatialMeasurements)
saveRDS(sce, "data/Head and Neck/sce.RDS")

