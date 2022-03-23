#' calculate_spatial_autocorrelation
#'
#' @description Calculate the spatial autocorrelation metric of the grid squares
#'   that have a certain pattern. -- I don't understand what a 'certain pattern' is, but maybe most people would find this obvious?
#'
#' @param raster_obj Raster object in the form of the output of raster function.
#' @param metric String. The method for calculating spatial autocorrelation.
#'   Choose from "globalmoran" and "GearyC". -- change code so that misspelling "GearyC" doesn't result in "globalmoran"
#' @export
#' @examples 
#' grid <- grid_metrics(SPIAT::defined_image, FUN = calculate_entropy, n_split = 5,
#' cell_types_of_interest=c("Tumour","Immune3"), feature_colname = "Cell.Type")
#' calculate_spatial_autocorrelation(grid, metric = "globalmoran")

calculate_spatial_autocorrelation <- function(raster_obj, metric = "globalmoran"){
  raster_obj@data@values[is.na(raster_obj@data@values)] = 0
  if (metric == "GearyC"){
    return(elsa::geary(raster_obj, d1=0, d2=600))
  }
  return(elsa::moran(raster_obj, d1=0, d2=600))
}
