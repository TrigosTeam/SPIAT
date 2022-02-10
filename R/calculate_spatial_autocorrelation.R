#' calculate_spatial_autocorrelaiton
#'
#' @description Calculate the spatial autocorrelation metric of the grid squares that have a certain pattern
#'
#' @param raster_obj Raster object in the form of the output of raster function
#' @param metric String The method for calculating spatial autocorrelation. Choose from "globalmoran" and "GearyC"
#' @import raster
#' @importFrom elsa geary moran
#' @export

calculate_spatial_autocorrelaiton <- function(raster, metric = "globalmoran"){
  raster@data@values[is.na(raster@data@values)] = 0
  if (metric == "GearyC"){
    return(geary(raster, d1=0, d2=600))
  }
  return(moran(raster, d1=0, d2=600))
}