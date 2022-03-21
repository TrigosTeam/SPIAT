#' compute_gradient
#'
#' @description Calculate the metrics for all of the specified radii	
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param radii Vector specifying the range of radii for the metrics to be calculated
#' @param FUN Variable name specifying the metric
#' @param ... Arguments of FUN
#' @return A list of the metrics under all radii
#' @export
#' 
#' @examples 
#' gradient_positions <- c(30, 50, 100)
#' gradient_entropy <- compute_gradient(SPIAT::defined_image, radii = gradient_positions, 
#' FUN = calculate_entropy,  cell_types_of_interest = c("Immune1","Immune2"),
#' feature_colname = "Cell.Type")

compute_gradient <- function(sce_object, radii, FUN, ...){
  list.metric <- list() 
  for (i in 1:length(radii)){
    metric <- FUN(sce_object,radius = radii[i], ...)
    list.metric[[i]] <- metric 
  }
  return(list.metric)
}
