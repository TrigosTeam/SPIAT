#' define_structure
#'
#' @description Identify the cells that compose the invasive front, the infiltration and the exclusion
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param names_of_immune_cells Vector indicating the names of potential immune cells
#' @param n_margin_layers Integer Specifying the number of layers of cells that compose the internal/external margins
#' @import SingleCellExperiment
#' @import dplyr
#' @export


define_structure <- function(sce_object, names_of_immune_cells, n_margin_layers = 5){
  
  # calculate the width of internal/external margin
  min_dist <- average_minimum_distance(sce_object)
  margin_dist <- n_margin_layers * min_dist
  print(paste("How many layers of cells in the external/internal margin:", n_margin_layers))
  print(paste("The width of internal/external margin:", margin_dist))
  
  #CHECK if the distance to bordering cells is calculated
  if (is.null(sce_object$Distance.To.Border)){
    stop(paste("Distance.To.Border not calculated yet for", deparse(substitute(sce_object))))
  }
  
  data <- data.frame(colData(sce_object))
  data[,"Structure"] <- data$Region
  data[intersect(which(data$Region == "Inside"),which(data$Cell.Type %in% names_of_immune_cells)), "Structure"] <- "Infiltrated.immune"
  data[intersect(which(data$Region == "Outside"),which(data$Cell.Type %in% names_of_immune_cells)), "Structure"] <- "Stromal.immune"
  data[intersect(which(data$Distance.To.Border < margin_dist), which(data$Region == "Inside")), "Structure"] <- "Internal.margin"
  data[intersect(which(data$Structure == "Internal.margin"), which(data$Cell.Type %in% names_of_immune_cells)), "Structure"] <- "Internal.margin.immune"
  data[intersect(which(data$Distance.To.Border < margin_dist), which(data$Region == "Outside")), "Structure"] <- "External.margin"
  data[intersect(which(data$Structure == "External.margin"), which(data$Cell.Type %in% names_of_immune_cells)), "Structure"] <- "External.margin.immune"
  
  colData(sce_object)$Structure <- data[,"Structure"]
  
  return(sce_object)
}
