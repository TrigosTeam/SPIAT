#' define_structure
#'
#' @description After identifying the bordering cells of tumour regions and
#'   calculating the distances of each cell to the tumour bordering cells, this
#'   function further identifies the cells that are located in the inside and
#'   outside of the tumour regions, and in the internal and external tumour
#'   margins. It also identifies the immune cells that are infiltrated, stromal,
#'   internal margin or external margin immune cells.
#'
#' @param sce_object SingleCellExperiment object whose metadata has the
#'   information of tumour bordering cells and cell distances to tumour border
#'   (has columns "Region" and "Distance.To.Border").
#' @param names_of_immune_cells String Vector of the names of immune cells.
#' @param n_margin_layers Integer. The number of layers of cells that compose
#'   the internal/external tumour margins.
#' @param feature_colname String Specifying the column that contains the names
#'   of the immune cells.
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
