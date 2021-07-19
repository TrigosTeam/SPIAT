#' average_nearest_neighbor_index
#'
#' @description Calculate the the ANN index of a specified type of cells. Return a list that contains 
#' the pattern of cells ("Clustered", "Dispersed", or "Random")  
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param reference_cell Cells positive for this marker will be used as reference
#' @param column String specify the selected column for reference_cell
#' @import SingleCellExperiment
#' @importFrom spatstat.geom nndist
#' @importFrom stats pnorm
#' @import ggplot2
#' @export


average_nearest_neighbor_index <- function(sce_object, reference_cell, column){
  
  ppp <- format_sce_to_ppp(sce_object)
  
  data <- data.frame(colData(sce_object))[,c(column ,"Cell.X.Position","Cell.Y.Position")]
  data <- data[which(data[,column] == reference_cell),c("Cell.X.Position","Cell.Y.Position") ]
  
  if(nrow(data) == 0){
    print("No reference cells found")
    ANN_index <- list(pattern=NA,`p-value`=NA)
  }else{
    object<-format_colData_to_sce(data)
    
    ppp_object <- format_sce_to_ppp(object)
    ann.p <- mean(nndist(ppp_object, k=1))
    
    n <- ppp_object$n      # Number of points
    
    x <- ppp$window$xrange[2] - ppp$window$xrange[1]
    y <- ppp$window$yrange[2] - ppp$window$yrange[1]
    area <- x*y
    ann.e <- 0.5/sqrt(n/area)
    se <- 0.26136/sqrt(n*n/area)
    z <- (ann.p - ann.e)/se
    p <- pnorm(-abs(z))
    
    if (p <= 5e-6){
      if (ann.p <= ann.e){
        pattern <- "Clustered"
      }
      else{pattern <- "Dispersed"}
    }
    else{pattern <- "Random"}
    
    ANN_index <- list(pattern, p)
    names(ANN_index) <- c("pattern" ,"p-value")
  }
  
  return(ANN_index)
}

