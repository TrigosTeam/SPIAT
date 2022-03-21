#' grid_metrics
#'
#' @description Calculate a specified metric for each little grid in the image and plot the metrics for the grids 
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param n_split Integer specifying the number of splits for the calculation of metrics
#' @param FUN Variable name specifying the metric to be calculated
#' @param ... arguments of FUN 
#' @return a list of the metrics of all grids
#' @export
#' @examples 
#' grid <- grid_metrics(SPIAT::defined_image, FUN = calculate_entropy, n_split = 5,
#' cell_types_of_interest=c("Tumour","Immune3"), feature_colname = "Cell.Type")

grid_metrics <- function(sce_object, FUN, n_split, ...){
  split <- image_splitter(sce_object,n_split)
  list.metric <- list()
  for (i in 1:length(split)){
    if(nrow(split[[i]]) > 0){
      sce <- try(quiet_basic(format_colData_to_sce(split[[i]])))
    }else{
      sce <- NULL
    }
    if (class(sce) == "SingleCellExperiment" || class(sce) == "SummarizedExperiment"){
      metric <-  quiet_basic(FUN(sce, ...))
      if (length(metric)==0){
        metric <- 0.0
      }
      list.metric[[i]] <- metric
    }
    else{
      list.metric[[i]] <- 0.0
    }
  }
  x <- raster::raster(ncol=n_split, nrow=n_split, xmn=0, xmx=max(sce_object$Cell.X.Position), ymn=0, 
              ymx=max(sce_object$Cell.Y.Position))
  raster::values(x) <- unlist(list.metric)
  y <- raster::flip(x, direction='y')
  raster::plot(y, main = paste("Plot ",as.character(substitute(FUN)), " of ", 
                       attr(sce_object, "name"), sep = ""))
  
  return(y)
}

