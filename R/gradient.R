#' gradient
#'
#' @description Calculate the metrics for all of the specified radii	
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param radii Vector specifying the range of radii for the metrics to be calculated
#' @param FUN Variable name specifying the metric
#' @param ... Arguments of FUN 
#' @importFrom SummarizedExperiment colData 
#' @return A list of the metrics under all radii
#' @export
#' 
gradient <- function(sce_object, radii, FUN, ...){
    
  list.metric <- list() 
  for (i in 1:length(radii)){
  
    metric <- FUN(sce_object,radius = radii[i], ...)
    list.metric[[i]] <- metric 
  }
  return(list.metric)
}
