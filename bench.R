library(profvis)
devtools::load_all()
library(tidySingleCellExperiment)

data(kerenSCE)


tumour <- c("Keratin_Tumour", "Tumour")

bcells <- c("B")
tcells <- c("CD3_Cell", "CD4_Cell", "CD8_Cell", "Tregs")
myeloid <- c("Dc/Mono", "DC", "Mono/Neu", "Macrophages", "Neutrophils")

endothelial <- c("Endothelial")
mesenchymal <- c("Mesenchymal")

tissue <- c(endothelial, mesenchymal)
immune <- c(bcells, myeloid, "NK", "other immune")

all <- c(tumour, tissue, immune, "Unidentified")

kerenSCE <- kerenSCE |>
    filter(imageID %in% c(5, 6)) |>
    filter(!(cellType %in% tcells))

parentDf <- parentCombinations(
    all = all,
    tumour,
    bcells,
    myeloid,
    endothelial,
    mesenchymal,
    tissue,
    immune
)


p <- profvis(
    kerenKontextual <- Kontextual(
        cells = kerenSCE,
        parentDf = parentDf,
        r = 100,
        cores = 1
    )
)

htmlwidgets::saveWidget(p, "profile.html")
