#' format_sce_to_ppp
#'
#' @description Formats an sce object into a ppp object
#' which has the x,y coordinates, phenotypes as markers 
#' window specifies the range of x and y coordinates
#'
#' @export
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @importFrom spatstat.geom ppp 
#' @importFrom grDevices chull

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


