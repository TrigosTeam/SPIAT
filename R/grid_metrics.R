#' grid_metrics
#'
#' @description Calculate a specified metric for each little grid in the image and plot the metrics for the grids 
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param n_split Integer specifying the number of splits for the calculation of metrics
#' @param FUN Variable name specifying the metric to be calculated
#' @param ... arguments of FUN 
#' @importFrom SummarizedExperiment colData 
#' @importFrom raster raster
#' @return a list of the metrics of all grids
#' @export

grid_metrics <- function(sce_object, FUN, n_split, ...){
  # split the image first
  split <- image_splitter(sce_object,n_split)
  
  # for each splitted image
  list.metric <- list()
  for (i in 1:length(split)){
    # make sure there are data in the splitted image
    try <- try(sce <- format_colData_to_sce(split[[i]]))
    if (class(try) == "SingleCellExperiment"){
      # make sure the calculation on the splitted image is feasible
      metric <- try(metric <- FUN(sce, ...), silent = T)
      # or the metric for this splitted image will be 0.0
      if (class(metric) == "try-error"){
        metric <- 0.0
      } 
      # mixing score (special case because mixing_score_multiple returns a list)
      else if (str_detect(deparse(substitute(FUN)), paste0("^", "mixing_score"))) {
        metric <- metric$Summary$Normalised_mixing_score
      }
    }
    # if there is no data in the splitted image, return 0.0
    else{
      metric <- 0.0
    }
    list.metric[[i]] <- metric
  }
  
  # create raster object
  x <- raster(ncol=n_split, nrow=n_split, xmn=0, xmx=max(sce_object$Cell.X.Position), ymn=0, 
              ymx=max(sce_object$Cell.Y.Position))
  # fill the metrics into the raster
  values(x) <- unlist(list.metric)
  # plot
  plot(x, main = paste("Plot", deparse(substitute(FUN)), "of", attr(sce_object, "name")))
  
  return(x)
}

