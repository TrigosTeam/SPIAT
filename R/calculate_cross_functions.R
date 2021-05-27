#' calculate_cross_functions
#'
#' @description calculate the cross functions between two cell types to see if there is mixing between them
#' 
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param method String that is the method for dependence calculation. Options: "Gcross", "Kcross", "Lcross", "Jcross"
#' @param phenotypes Vector of phenotypes of interest
#' @param column String that is the name of the column of the types
#' @import SingleCellExperiment
#' @importFrom spatstat.core Gcross Kcross Lcross Jcross
#' @importFrom spatstat.geom ppp 
#' @export 

calculate_cross_functions <- function(sce_object, method = "Kcross", 
                                 phenotypes, column) {
  
  formatted_data <- colData(sce_object)
  
  #CHECK
  if (!all(phenotypes %in% formatted_data[[column]])) {
    stop("Cell type not found!")
  }
  
  # get x, y coordinates and specified cell types from formatted_data
  x <- formatted_data$Cell.X.Position
  y <- formatted_data$Cell.Y.Position
  marks <- formatted_data[[column]]
  
  # get windows
  x_summary <- summary(formatted_data$Cell.X.Position)
  x_min <- as.numeric(x_summary[1])
  x_max <- as.numeric(x_summary[6])
  y_summary <- summary(formatted_data$Cell.Y.Position)
  y_min <- as.numeric(y_summary[1])
  y_max <- as.numeric(y_summary[6])
  
  
  # format sce to ppp
  ppp_object <- ppp(x, y, c(x_min,x_max),c(y_min,y_max), marks = marks)
  ppp_object$marks <- as.factor(ppp_object$marks)
  
  if (method == "Gcross"){
    p <- Gcross(ppp_object, phenotypes[1],phenotypes[2],correction = "border")
    plot(p, main = paste("cross G function",attr(sce_object,"name")))
  }
  if (method == "Kcross"){
    p <- Kcross(ppp_object, phenotypes[1],phenotypes[2],correction = "border")
    plot(p, main = paste("cross K function",attr(sce_object,"name")))
  }
  if (method == "Lcross"){
    p <- Lcross(ppp_object, phenotypes[1],phenotypes[2],correction = "border")
    plot(p, main = paste("cross L function",attr(sce_object,"name")))
  }
  if (method == "Jcross"){
    p <- Jcross(ppp_object, phenotypes[1],phenotypes[2],correction = "border")
    plot(p, main = paste("cross J function",attr(sce_object,"name")))
  }
  return(p)
}

