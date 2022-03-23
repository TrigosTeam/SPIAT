#' define_structure
#'
#' @description Identify the cells that compose the invasive front, the infiltration and the exclusion -- elaborate on what these mean? Or cite your paper for definitions?
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param names_of_immune_cells Vector indicating the names of potential immune cells
#' @param n_margin_layers Integer Specifying the number of layers of cells that compose the internal/external margins
#' @param feature_colname String Specifying which column the names of immune cells are under
#' @import dplyr
#' @export
#' @examples 
#' sce_border <- identify_bordering_cells(SPIAT::defined_image, 
#' reference_cell = "Tumour", feature_colname = "Cell.Type", n_to_exclude = 10)
#' sce_dist <- calculate_distance_to_tumour_margin(sce_border)
#' sce_structure <- define_structure(sce_dist, 
#' names_of_immune_cells = c("Immune1","Immune2","Immune3"),
#' feature_colname = "Cell.Type", n_margin_layers = 5)
#' plot_cell_categories(sce_structure, feature_colname = "Structure")

define_structure <- function(sce_object, names_of_immune_cells, 
                             feature_colname = "Cell.Type",
                             n_margin_layers = 5){
  
  # calculate the width of internal/external margin
  min_dist <- average_minimum_distance(sce_object)
  margin_dist <- n_margin_layers * min_dist
  print(paste("How many layers of cells in the external/internal margin:", n_margin_layers))
  print(paste("The width of internal/external margin:", margin_dist))
  
  #CHECK if the distance to bordering cells is calculated
  if (is.null(sce_object$Distance.To.Border)){
    stop(paste("Distance.To.Border not calculated yet for", deparse(substitute(sce_object))))
  }
  
  data <- data.frame(SummarizedExperiment::colData(sce_object))
  data[,"Structure"] <- data$Region
  data[intersect(which(data$Region == "Inside"),which(data[[feature_colname]] %in% names_of_immune_cells)), "Structure"] <- "Infiltrated.immune"
  data[intersect(which(data$Region == "Outside"),which(data[[feature_colname]] %in% names_of_immune_cells)), "Structure"] <- "Stromal.immune"
  data[intersect(which(data$Distance.To.Border < margin_dist), which(data$Region == "Inside")), "Structure"] <- "Internal.margin"
  data[intersect(which(data$Structure == "Internal.margin"), which(data[[feature_colname]] %in% names_of_immune_cells)), "Structure"] <- "Internal.margin.immune"
  data[intersect(which(data$Distance.To.Border < margin_dist), which(data$Region == "Outside")), "Structure"] <- "External.margin"
  data[intersect(which(data$Structure == "External.margin"), which(data[[feature_colname]] %in% names_of_immune_cells)), "Structure"] <- "External.margin.immune"
  
  SummarizedExperiment::colData(sce_object)$Structure <- data[,"Structure"]
  
  return(sce_object)
}
