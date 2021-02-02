#' calculate_k_function 
#' 
#' @description Returns K function of cell distribution (inhomogenous distributions)
#' 
#' @param sce_object SingleCellExperiment object used as input data for K function 
#' @param determine_confidence indicate the confidence interval of the K function calculated
#' 

calculate_k_function <- function(sce_object, cell_phenotypes_of_interest = NULL, determine_confidence = NULL, plot_envelope = NULL){
  #import spatstat package
  library(spatstat)
  
  #create point pattern object
  point_pattern <- ppp(colData(sce_object)[,2], colData(sce_object)[,3], c(min(colData(sce_object)[,2]), max(colData(sce_object)[,2])), c(min(colData(sce_object)[,3]), max(colData(sce_object)[,3])), marks = as.factor(colData(sce_object)[,1]))
  
  if(!is.null(cell_phenotypes_of_interest)) {
    subset_of_interest <- subset(point_pattern, marks == cell_phenotypes_of_interest)
    point_pattern <- subset_of_interest
  } 
  k_function <- Kinhom(point_pattern)
  

  
  if (!is.null(determine_confidence)) {
    confidence <- lohboot(point_pattern, "Kinhom")
    plot(confidence)
  }
  if (!is.null(plot_envelope)) {
    envelope_k <- envelope(point_pattern, Kinhom, nsim= 49)
    plot(envelope_k, main = "K-function")
  } 
  else {
    plot(k_function, main = "K-function")
  }
  return(k_function)
}