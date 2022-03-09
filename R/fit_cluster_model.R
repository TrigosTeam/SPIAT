#' Fit Cauchy model to point pattern
#' 
#' fit_cluster_model fits a Cauchy cluster model to the data given in point_pattern. 
#' If phenotypes_of_interest is given only a subset of data will be used to create the model. 
#' 
#' Inputs - 
#' point_pattern: ppp object representing data can be marked or unmarked
#' phenotypes_of_interest: String with names of phenotypes of interest, if specified the model will be produced based off the subset data. Set to NULL by default
#' 
#' Output - 
#' model_parameters: List containing model parameters, kappa, mean cluster size and estimated intensity
#' 
#' @param point_pattern ppp object representing data, can be marked point pattern of whole sample or of a single phenotype
#' @param phenotypes_of_interest specifies which phenotype to base cluster model
#' 
#' @importFrom spatstat.core kppm 
#' @importFrom spatstat.geom unmark

fit_cluster_model <- function(point_pattern, phenotypes_of_interest = NULL){
  # Check for specified cell phenotype and subset if necessary 
  if(!is.null(phenotypes_of_interest)) {
    temp <- find_subset(point_pattern, cell_phenotypes_of_interest = phenotypes_of_interest)
    point_pattern <- temp
  } 
  # Unmark point_pattern 
  point_pattern <- unmark(point_pattern)
  # Fit Cauchy cluster model
  model <- kppm(point_pattern, clusters = "Cauchy")
  #Extract and return relevant data
  model_parameters <- list(kappa = model$clustpar[1], mean_cluster_size = model$mu, estimated_intensity = model$lambda)
   
  return(model_parameters)
}
