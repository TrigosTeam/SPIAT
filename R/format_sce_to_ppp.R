#' Create point pattern object from single cell experiment object
#' 
#' Input - Single cell experiment object to be transformed into ppp object,
#'  marker phenotypes will be allocated as marks to ppp object
#' 
#' Output - ppp object, marked to reflect initial image
#' 
#' @param sce_object sce object to be transformed into ppp object


format_sce_to_ppp <- function(sce_object){
  # insert sce_object data into ppp object
  point_pattern <- ppp(colData(sce_object)[,2], colData(sce_object)[,3], c(min(colData(sce_object)[,2]), max(colData(sce_object)[,2])), c(min(colData(sce_object)[,3]), max(colData(sce_object)[,3])), marks = as.factor(colData(sce_object)[,1]))
  
  #return ppp object
  return(point_pattern)
}
=======
#' format_sce_to_ppp
#'
#' @description Formats an sce object into a ppp object
#' which has the x,y coordinates, phenotypes as markers 
#' window specifies the range of x and y coordinates
#'
#' @export
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce

format_sce_to_ppp <- function(sce_object) {
  
  # get x, y coordinates and phenotypes from sce object
  sce_data <- colData(sce_object)
  x <- sce_data$Cell.X.Position
  y <- sce_data$Cell.Y.Position
  marks <- sce_data$Phenotype
  
  # get windows
  x_summary <- summary(sce_data$Cell.X.Position)
  x_min <- as.numeric(x_summary[1])
  x_max <- as.numeric(x_summary[6])
  y_summary <- summary(sce_data$Cell.Y.Position)
  y_min <- as.numeric(y_summary[1])
  y_max <- as.numeric(y_summary[6])
  
  # get ploy window
  X <- data.frame(x,y)
  hpts <- chull(X)
  poly_window <- list(x=rev(X[hpts, 1]), y=rev(X[hpts, 2]))
  
  # format sce to ppp
  ppp_object <- ppp(x, y, poly = poly_window, 
                    marks = marks)
  
  return(ppp_object)
}


