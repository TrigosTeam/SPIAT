#' define_structure
#'
#' @description Identify the cells that compose the invasive front, the infiltration and the exclusion
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param names_of_immune_cells vector indicating the names of potential immune cells
#' @param n_invasive Number specifying the number of cells that are included in the invasive front
#' @import SingleCellExperiment
#' @import dplyr
#' @export


define_structure <- function(sce_object, names_of_immune_cells, n_invasive = 5){
  
  # find the distance of invasive front
  min_dist <- average_minimum_distance(sce_object)
  invasive_dist <- n_invasive * min_dist
  print(invasive_dist)
  
  #CHECK if the distance to bordering cells is calculated
  if (is.null(sce_object$Distance.To.Border)){
    stop(paste("Distance.To.Border not calculated yet for", deparse(substitute(sce_object))))
  }
  
  data <- data.frame(colData(sce_object))
  data[,"Structure"] <- data$Region
  data[intersect(which(data$Region == "In"),which(data$Cell.Type %in% names_of_immune_cells)), "Structure"] <- "Infiltration"
  data[intersect(which(data$Region == "Out"),which(data$Cell.Type %in% names_of_immune_cells)), "Structure"] <- "Exclusion"
  data[intersect(which(data$Distance.To.Border < invasive_dist), which(data$Region == "In")), "Structure"] <- "Invasive.front.in"
  data[intersect(which(data$Structure == "Invasive.front.in"), which(data$Cell.Type %in% names_of_immune_cells)), "Structure"] <- "I.f.immune.in"
  data[intersect(which(data$Distance.To.Border < invasive_dist), which(data$Region == "Out")), "Structure"] <- "Invasive.front.out"
  data[intersect(which(data$Structure == "Invasive.front.out"), which(data$Cell.Type %in% names_of_immune_cells)), "Structure"] <- "I.f.immune.out"
  
  colData(sce_object)$Structure <- data[,"Structure"]
  
  return(sce_object)
}
