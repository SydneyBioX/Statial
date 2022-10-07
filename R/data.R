#' Head and neck cutaneous squamous cell carcinoma IMC data 
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
#' This is a single image obtained from the simulateCompartment function obtained using the following code:
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


#' Patient 6 image from Keren et al
#'
#'
#' This is a single image obtained from the simulateCompartment function obtained using the following code:
#' 
#'
#' @format headData a SingleCellExperiment object
#' @aliases 
#' kerenImage
"kerenImage"