#' Split an image into grid and calculates a metric for each grid square
#'
#' @description Calculates a specified metric for each grid tile in the image
#'   and plots the metrics for the grid tiles.
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param n_split Integer specifying the number of splits for the calculation of
#'   metrics. This number is the splits on each side (e.g. `n_split` = 3 means
#'   the image will be split into 9 tiles.)
#' @param FUN Variable name specifying the metric to be calculated.
#' @param ... Arguments of FUN
#' @return A list of the metrics of all grid tiles
#' @export
#' @examples
#' grid <- grid_metrics(SPIAT::defined_image, FUN = calculate_entropy, n_split = 5,
#' cell_types_of_interest=c("Tumour","Immune3"), feature_colname = "Cell.Type")

grid_metrics <- function(spe_object, FUN, n_split, ...){
  split <- image_splitter(spe_object,n_split)
  list.metric <- list()
  for (i in seq_len(length(split))){
    spe <- split[[i]]
    if(nrow(spe) == 0){spe <- NULL}
    if (methods::is(spe,"SpatialExperiment")){
      metric <-  quiet_basic(FUN(spe, ...))
      if (length(metric)==0){
        metric <- 0.0
      }
      list.metric[[i]] <- metric
    }
    else{
      list.metric[[i]] <- 0.0
    }
  }
  x <- raster::raster(ncol=n_split, nrow=n_split, xmn=0, ymn=0, 
                      xmx=max(SpatialExperiment::spatialCoords(spe_object)[,"Cell.X.Position"]), 
                      ymx=max(SpatialExperiment::spatialCoords(spe_object)[,"Cell.X.Position"]))
  raster::values(x) <- unlist(list.metric)
  y <- raster::flip(x, direction='y')
  raster::plot(y, main = paste("Plot ",attr(spe_object, "name"), " ",as.character(substitute(FUN)), 
                       sep = ""))
  
  return(y)
}

