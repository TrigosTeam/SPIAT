#' define_regions
#'
#' @description Identify the cells that compose the invasive front, the infiltration and the exclusion
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param names_of_immune_cells vector indicating the names of potential immune cells
#' @param n_invasive Number specifying the number of cells that are included in the invasive front
#' @import SingleCellExperiment
#' @import dplyr
#' @export


define_regions <- function(sce_object, names_of_immune_cells = c("CD4","CD8"), n_invasive = 5){
  
  # find the distance of invasive front
  min_dist <- average_minimum_distance(sce_object)
  invasive_dist <- n_invasive * min_dist
  print(invasive_dist)
  
  #CHECK if the distance to bordering cells is calculated
  if (is.null(sce_object$Distance.To.Border)){
    stop(paste("Distance.To.Border not calculated yet for", deparse(substitute(sce_object))))
  }
  
  data <- data.frame(colData(sce_object))
  data[intersect(which(data$Region == "In"),which(data$Cell.Type %in% names_of_immune_cells)), "Region2"] <- "Infiltration"
  data[intersect(which(data$Region == "Out"),which(data$Cell.Type %in% names_of_immune_cells)), "Region2"] <- "Exclusion"
  data[intersect(which(data$Distance.To.Border < invasive_dist), which(data$Region == "In")), "Region2"] <- "Invasive.front.in"
  data[intersect(which(data$Distance.To.Border < invasive_dist), which(data$Region == "Out")), "Region2"] <- "Invasive.front.out"
  
  colData(sce_object)$Region2 <- data[,"Region2"]
  
  return(sce_object)
}
