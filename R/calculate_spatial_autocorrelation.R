#' calculate_spatial_autocorrelation
#'
#' @description Calculate the spatial autocorrelation metric of the grid squares
#'   that have a certain pattern.
#'
#' @param raster_obj Raster object in the form of the output of raster function.
#' @param metric String. The method for calculating spatial autocorrelation.
#'   Choose from "globalmoran" and "GearyC".
#' @export

calculate_spatial_autocorrelaiton <- function(raster_obj, metric = "globalmoran"){
  raster_obj@data@values[is.na(raster_obj@data@values)] = 0
  if (metric == "GearyC"){
    return(elsa::geary(raster_obj, d1=0, d2=600))
  }
  return(elsa::moran(raster_obj, d1=0, d2=600))
}
