#' calculate_spatial_autocorrelation
#'
#' @description Takes the result of \code{\link{grid_metrics}} (a RasterLayer
#'   object) and calculate its spatial autocorrelation.
#'
#' @param raster_obj Raster object in the form of the output of
#'   \code{\link{grid_metrics}}.
#' @param metric String. The method for calculating spatial autocorrelation.
#'   Choose from "globalmoran" and "GearyC".
#' @export
#' @return A number is returned
#' @examples
#' grid <- grid_metrics(SPIAT::defined_image, FUN = calculate_entropy, 
#' n_split = 5, cell_types_of_interest=c("Tumour","Immune3"), 
#' feature_colname = "Cell.Type")
#' calculate_spatial_autocorrelation(grid, metric = "globalmoran")

calculate_spatial_autocorrelation <- function(raster_obj, 
                                              metric = "globalmoran"){
    raster_obj@data@values[is.na(raster_obj@data@values)] = 0
    if (metric == "GearyC"){
        return(elsa::geary(raster_obj, d1=0, d2=600))
    }
    return(elsa::moran(raster_obj, d1=0, d2=600))
}
