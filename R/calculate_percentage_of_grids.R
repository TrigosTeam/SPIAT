#' calculate_percentage_of_grids
#'
#' @description Calculate the percentage of the grid squares that have a certain
#'   pattern. --might need to contextualise by describing how this is used with other functions?
#'
#' @param raster_obj Raster object in the form of the output of `grid_metrics`
#'   function.
#' @param threshold Numeric. The threshold for defining the pattern.
#' @param above Boolean. Indicating whether the pattern is above (TRUE) or below
#'   (FALSE) the threshold.
#' @export
#' @examples
#' grid <- grid_metrics(SPIAT::defined_image, FUN = calculate_entropy, n_split = 5,
#' cell_types_of_interest=c("Tumour","Immune3"), feature_colname = "Cell.Type")
#' calculate_percentage_of_grids(grid, threshold = 0.75, above = TRUE)

calculate_percentage_of_grids <- function(raster_obj, threshold, above){
  raster_obj@data@values[is.na(raster_obj@data@values)] <- 0
  if (above) p <- sum(raster_obj@data@values>=threshold)/length(raster_obj@data@values) * 100
  else p <- sum(raster_obj@data@values<threshold)/length(raster_obj@data@values) * 100
  return(p)
}
