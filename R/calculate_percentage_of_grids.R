#' calculate_percentage_of_grids
#'
#' @description Calculate the percentage of the grid squares that have a certain pattern
#'
#' @param raster_obj Raster object in the form of the output of raster function
#' @param threshold Number that is the threshold for defining the pattern
#' @param above Boolean idicating the pattern is above (TRUE) or below (FALSE) the threshold
#' @export

calculate_percentage_of_grids <- function(raster_obj, threshold, above){
  raster_obj@data@values[is.na(raster_obj@data@values)] <- 0
  if (above) p <- sum(raster_obj@data@values>=threshold)/length(raster_obj@data@values) * 100
  else p <- sum(raster_obj@data@values<threshold)/length(raster_obj@data@values) * 100
  return(p)
}