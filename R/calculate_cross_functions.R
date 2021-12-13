#' calculate_cross_functions
#'
#' @description calculate the cross functions between two specified cell types 
#' 
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param method String that is the method for dependence calculation. Options: "Gcross", "Kcross", "Lcross", "Jcross"
#' @param phenotypes Vector of phenotypes of interest
#' @param feature_colname String that is the name of the column of the types
#' @param plot_results TRUE if result to be plotted, FALSE if not. In either case, an object with the results is returned
#' @importFrom spatstat.core Gcross Kcross Lcross Jcross
#' @importFrom spatstat.geom ppp 
#' @export 

calculate_cross_functions <- function(sce_object, method = "Kcross", 
                                      phenotypes, feature_colname, plot_results = T) {
  
  #CHECK
  if (!all(phenotypes %in% formatted_data[[feature_colname]])) {
    stop("Cell type not found!")
  }
  
  # format sce to ppp object
  ppp_object <- format_sce_to_ppp(sce_object)
  ppp_object$marks <- as.factor(ppp_object$marks)
  
  if (method == "Gcross"){
    p <- Gcross(ppp_object, phenotypes[1],phenotypes[2],correction = "border")
    if(plot_results){
      plot(p, main = paste("cross G function",attr(sce_object,"name")))
    }
  }
  if (method == "Kcross"){
    p <- Kcross(ppp_object, phenotypes[1],phenotypes[2],correction = "border")
    if(plot_results){
      plot(p, main = paste("cross K function",attr(sce_object,"name")))
    }
  }
  if (method == "Lcross"){
    p <- Lcross(ppp_object, phenotypes[1],phenotypes[2],correction = "border")
    if(plot_results){
      plot(p, main = paste("cross L function",attr(sce_object,"name")))
    }
  }
  if (method == "Jcross"){
    p <- Jcross(ppp_object, phenotypes[1],phenotypes[2],correction = "border")
    if(plot_results){
      plot(p, main = paste("cross J function",attr(sce_object,"name")))
    }
  }
  return(p)
}

