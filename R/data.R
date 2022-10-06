#' Head and Neck cancer IMC data 
#'
#'
#' This is a imaging mass cytometry dataset from Ferguson et al 2022. The object contains 29 images.
#'
#'
#' @format headData a SingleCellExperiment object
#' @aliases 
#' headSCE
"headSCE"



#' Example simulated image
#'
#'
#' This is a single image obtained from the simulateCompartment function using the following code:
#' 
#' set.seed(10)
#' 
#' #Simulating example image
#' simulation = simulateCompartment(includeTissue = FALSE)
#' 
#' #Selecting image where a significant conditional relationship exists
#' exampleImage = simulation$sig
#'
#' @format headData a SingleCellExperiment object
#' @aliases 
#' exampleImage
"exampleImage"


