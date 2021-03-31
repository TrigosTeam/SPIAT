#' number_of_cells_within_radius
#'
#' @description Calculates the number of cells positive for a target marker within a pre-defined radius around cells positive for a reference marker

#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param reference_marker Marker to be used for reference cells
#' @param target_marker Marker to be used for target cells
#' @param radius Radius around the reference cells
#' @importFrom dbscan frNN
#' @return The median number of target cells within the specified radius

number_of_cells_within_radius <- function (sce_object, reference_marker, target_marker, radius = 20) 
{
  formatted_data <- data.frame(colData(sce_object))
  formatted_data <- formatted_data %>% rownames_to_column("Cell.ID")
  expression_matrix <- assay(sce_object)
  markers <- rownames(expression_matrix)
  cell_ids <- colnames(expression_matrix)
  rownames(expression_matrix) <- NULL
  colnames(expression_matrix) <- NULL
  expression_matrix_t <- t(expression_matrix)
  expression_df <- data.frame(expression_matrix_t)
  colnames(expression_df) <- markers
  formatted_data <- cbind(formatted_data, expression_df)
  formatted_data <- formatted_data[complete.cases(formatted_data),]
  reference_cells <- formatted_data[grepl(reference_marker, 
                                          formatted_data$Phenotype), ]
  if (nrow(reference_cells) == 0) {
    return(NA)
  }
  target_cells <- formatted_data[grepl(target_marker, formatted_data$Phenotype),]
  if (nrow(target_cells) == 0) {
    return(NA)
  }

  common_cells <- reference_cells$Cell.ID[reference_cells$Cell.ID %in% 
                                            target_cells$Cell.ID]
  reference_cells <- reference_cells[!(reference_cells$Cell.ID %in% 
                                         common_cells), ]
  target_cells <- target_cells[!(target_cells$Cell.ID %in% 
                                   common_cells), ]
  reference_cell_cords <- reference_cells[, c("Cell.X.Position", 
                                              "Cell.Y.Position")]
  target_cell_cords <- target_cells[, c("Cell.X.Position", 
                                        "Cell.Y.Position")]
  search_result <- frNN(target_cell_cords, eps = radius, query = reference_cell_cords, 
                        sort = FALSE)
  number_of_target_by_reference <- sapply(search_result$id, length)

  return(median(number_of_target_by_reference[number_of_target_by_reference != 0]))
}
