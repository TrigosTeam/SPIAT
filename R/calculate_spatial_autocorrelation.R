#' calculate_spatial_autocorrelation
#'
#' @description Takes the result of \code{\link{grid_metrics}} (a RasterLayer
#'   object) and calculate its spatial autocorrelation.
#'
#' @param raster_obj Raster object in the form of the output of
#'   \code{\link{grid_metrics}}.
#' @param metric String. The method for calculating spatial autocorrelation.
#'   Choose from "globalmoran" and "GearyC".
#' @param d Numeric. Upper bound local distance. The argument `d2` from function 
#'   \link[elsa]{moran}. Default is NULL and the distance will be calculated 
#'   automatically from the number of splits and the extent of the grid image. 
#' @export
#' @return A number is returned
#' @examples
#' grid <- grid_metrics(SPIAT::defined_image, FUN = calculate_entropy, 
#' n_split = 5, cell_types_of_interest=c("Tumour","Immune3"), 
#' feature_colname = "Cell.Type")
#' calculate_spatial_autocorrelation(grid, metric = "globalmoran")

calculate_spatial_autocorrelation <- function(raster_obj, 
                                              metric = "globalmoran",
                                              d = NULL){
    requireNamespace("elsa", quietly = TRUE)
    extent <- min(raster_obj@extent@ymax, raster_obj@extent@xmax)
    if (is.null(d)) d <- round(extent/raster_obj@ncols * 2)
    raster_obj@data@values[is.na(raster_obj@data@values)] <- 0
    if (metric == "GearyC"){
        return(elsa::geary(raster_obj, d1=0, d2=d))
    }
    return(elsa::moran(raster_obj, d1=0, d2=d))
}
