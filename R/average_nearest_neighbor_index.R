#' ANN index for point pattern (clustering or dispersion)
#'
#' @description Calculate the the ANN index of a specified type of cells. The
#'   index indicates the clustering effect of a point pattern. It can be
#'   clustering, random or dispersion.
#'
#' @param sce_object SingleCellExperiment object in the form of the output of
#'   format_image_to_sce.
#' @param reference_cell String. Cells of this type are used as reference cells.
#' @param feature_colname String. Specify the selected column for
#'   `reference_cell`.
#' @export


average_nearest_neighbor_index <- function(sce_object, reference_cell, feature_colname){
  
  ppp <- format_sce_to_ppp(sce_object)
  formatted_data <- data.frame(SummarizedExperiment::colData(sce_object))
  
  data <- formatted_data[,c(feature_colname,"Cell.X.Position","Cell.Y.Position")]
  data <- data[which(data[,feature_colname] == reference_cell),
               c("Cell.X.Position","Cell.Y.Position") ]
  
  if(nrow(data) == 0){
    print("No reference cells found")
    ANN_index <- list(pattern=NA,`p-value`=NA)
  }else{
    object<-format_colData_to_sce(data)
    
    ppp_object <- format_sce_to_ppp(object)
    ann.p <- mean(spatstat.geom::nndist(ppp_object, k=1))
    
    n <- ppp_object$n # Number of points
    
    x <- ppp$window$xrange[2] - ppp$window$xrange[1]
    y <- ppp$window$yrange[2] - ppp$window$yrange[1]
    area <- x*y
    ann.e <- 0.5/sqrt(n/area)
    se <- 0.26136/sqrt(n*n/area)
    z <- (ann.p - ann.e)/se
    p <- stats::pnorm(-abs(z))
    
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

