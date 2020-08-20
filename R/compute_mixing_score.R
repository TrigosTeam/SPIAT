#' compute_mixing_score
#'
#' @description Returns the mixing score between a reference marker and a target
#' marker. The mixing score is defined by:
#' the number of target-reference interactions/number of reference-reference interactions
#' within a specified radius.
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param reference_marker String specifying the reference marker
#' @param target_marker String specifying the target marker
#' @param radius Integer specifying the radius. Only cells within this radius will be considered.
#' @import dplyr
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble rownames_to_column
#' @importFrom stats complete.cases
#' @importFrom dbscan frNN
#' @export

compute_mixing_score <- function(sce_object, reference_marker, target_marker, radius = 20) {

  formatted_data <- data.frame(colData(sce_object))
  formatted_data <- formatted_data[complete.cases(formatted_data),]
  formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

  #select the cells that contain the target marker
  target_cells <- formatted_data[grepl(target_marker, formatted_data$Phenotype), ]

  #select the cells that contain the reference marker
  reference_cells <- formatted_data[grepl(reference_marker, formatted_data$Phenotype), ]

  #Exclude cells that express both markers
  target_cells <- target_cells[!grepl(reference_marker, target_cells$Phenotype), ]
  reference_cells <- reference_cells[!grepl(target_marker, reference_cells$Phenotype), ]

  if (nrow(reference_cells) == 0) {
    stop("There are no reference cells of specified marker, cannot calculate mixing score")
  }
  if (nrow(target_cells) == 0) {
    stop("There are no target cells of specified marker")
  }

  print(paste("Number of reference cells:", nrow(reference_cells)))
  print(paste("Number of target cells:", nrow(target_cells)))

  reference_cell_cords <- reference_cells[,c("Cell.X.Position", "Cell.Y.Position")]
  target_cell_cords <- target_cells [,c("Cell.X.Position", "Cell.Y.Position")]

  #calculate the number of interactions between reference and target cells
  reference_target_result <- frNN(target_cell_cords, eps = radius, query = reference_cell_cords, sort = FALSE)
  #sum up the number of interactions
  reference_target_interactions <- sum(rapply(reference_target_result$id, length))
  print(paste("Number of reference-target interactions:", reference_target_interactions))

  #calculate the number of interactions between reference cells
  reference_reference_result <- frNN(reference_cell_cords, eps=radius, sort = FALSE)
  #sum up the number of interactions
  reference_reference_interactions <- sum(rapply(reference_reference_result$id, length))
  print(paste("Number of reference-reference interactions:", reference_reference_interactions))
  if (reference_reference_interactions != 0) {
    mixing_score <- reference_target_interactions/reference_reference_interactions
    return (mixing_score)
  } else {
    stop("There are no interactions between reference cells for the specified radius, cannot calculate mixing score")
  }
}
