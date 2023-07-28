#' Head and neck cutaneous squamous cell carcinoma IMC data
#'
#' @description
#' This is a subset of imaging mass cytometry dataset of head and neck cutaneous
#' squamous cell carcinoma from Ferguson et al 2022. The object contains 5 images.
#'
#' @usage data(headSCE)
#'
#' @format headData a SingleCellExperiment object
#' @aliases
#' headSCE
#'
#' @references Ferguson, A. L., Sharman, A. R., Allen, R. O., Ye, T., Lee, J. H.,
#' Low, T.-H. H., Ch'ng, S., Palme, C. E., Ashford, B., Ranson, M.,
#' Clark, J. R., Patrick, E., Gupta, R., & Palendira, U. (2022).
#' High-Dimensional and Spatial Analysis Reveals Immune Landscape–Dependent
#' Progression in Cutaneous Squamous Cell Carcinoma. Clinical Cancer Research,
#' OF1-OF12. ([DOI](https://doi.org/10.1158/1078-0432.Ccr-22-1332))
#'
"headSCE"



#' MIBI-TOF Breast cancer image
#'
#' @description
#' This is a single MIBI-TOF image of breast cancer from patient 6 of the
#' Keren et al 2018 dataset.
#'
#'
#' @usage data(kerenImage)
#'
#' @format headData a SingleCellExperiment object
#' @aliases
#' kerenImage
#'
#' @references Keren, L., Bosse, M., Marquez, D., Angoshtari, R., Jain,
#' S., Varma, S., Yang, S. R., Kurian, A., Van Valen, D., West, R., Bendall,
#' S. C., & Angelo, M. (2018). A Structured Tumor-Immune Microenvironment in
#' Triple Negative Breast Cancer Revealed by Multiplexed Ion Beam Imaging.
#' Cell, 174(6), 1373-1387.e1319. ([DOI](https://doi.org/10.1016/j.cell.2018.08.039))
#'
#'
"kerenImage"


#' MIBI-TOF Breast cancer intensities
#'
#' @description
#' This is a single MIBI-TOF data of breast cancer from patient 6 of the
#' Keren et al 2018 dataset.
#'
#'
#' @usage data(kerenSCE)
#'
#' @format headData a SingleCellExperiment object
#' @aliases
#' kerenSCE
#'
#' @references Keren, L., Bosse, M., Marquez, D., Angoshtari, R., Jain,
#' S., Varma, S., Yang, S. R., Kurian, A., Van Valen, D., West, R., Bendall,
#' S. C., & Angelo, M. (2018). A Structured Tumor-Immune Microenvironment in
#' Triple Negative Breast Cancer Revealed by Multiplexed Ion Beam Imaging.
#' Cell, 174(6), 1373-1387.e1319. ([DOI](https://doi.org/10.1016/j.cell.2018.08.039))
#'
#'
"kerenSCE"

#' Results from getStateChanges with the mixed effects model for kerenSCE
#'
#' Results from the call:
#' SCE <- getDistances(kerenSCE, nCores = 20)
#' SCE <- getAbundances(SCE, nCores = 20)
#' 
#' stateChangesMixed <- getStateChanges(
#'   singleCellData = SCE,
#'   Rs = c(200),
#'   typeAll = c("dist200", "abundance200"),
#'   method = "lm",
#'   isMixed = FALSE,
#'   nCores = 40)
#'
#' @format stateChanges a dataframe
#' @aliases 
#' stateChanges
"stateChangesFiltered"



#' Results from getStateChanges with the mixed effects model for kerenSCE
#'
#' Results from the call:
#' SCE <- getDistances(kerenSCE, nCores = 20)
#' SCE <- getAbundances(SCE, nCores = 20)
#' 
#' stateChangesMixed <- getStateChanges(
#'   singleCellData = SCE,
#'   Rs = c(200),
#'   typeAll = c("dist200", "abundance200"),
#'   method = "lm",
#'   isMixed = TRUE,
#'   nCores = 40)
#'
#' @format stateChangesMixed a dataframe
#' @aliases 
#' stateChangesMixed
"stateChangesMixed"


#' Results from getStateChanges with the mixed effects model when using
#' contamination as a covariate for kerenSCE.
#'
#' Results from the call:
#' SCE <- getDistances(kerenSCE, nCores = 20)
#' SCE <- getAbundances(SCE, nCores = 20)
#' SCE <- calcContamination(SCE)
#' 
#' stateChangesContam <- getStateChanges(
#'   singleCellData = SCE,
#'   Rs = c(200),
#'   typeAll = c("dist200", "abundance200"),
#'   covariates = c("rfMainCellProb"),
#'   method = "lm",
#'   isMixed = TRUE,
#'   nCores = 40)
#'
#' @format stateChangesContam a dataframe
#' @aliases 
#' stateChangesContam
"stateChangesContam"

