#' number_of_cells_within_radius
#'
#' @description Calculates the number of cells positive for a target marker within a pre-defined radius around cells positive for a reference marker

#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param reference_marker Marker to be used for reference cells
#' @param target_marker Marker to be used for target cells
#' @param radius Radius around the reference cells
#' @param column String specifying the column with the desired cell type annotations
#' @importFrom dbscan frNN
#' @return List of dataframes with the number of target cells of each of the reference cells
#' @export

number_of_cells_within_radius <- function (sce_object, reference_marker, target_marker, radius = 20,column) 
{
  formatted_data <- data.frame(colData(sce_object))
  formatted_data <- formatted_data %>% rownames_to_column("Cell.ID")
  formatted_data <- formatted_data[complete.cases(formatted_data),]
  all.df <- list()
  for (i in reference_marker) {
    reference_cells <- formatted_data[which(formatted_data[,column] == i), ]
    reference_cell_cords <- reference_cells[, c("Cell.ID", "Cell.X.Position","Cell.Y.Position")]
    dataframe <- remove_rownames(reference_cell_cords)
    dataframe <- dataframe %>% column_to_rownames("Cell.ID")
    reference_cell_cords <- reference_cells[, c( "Cell.X.Position","Cell.Y.Position")]
    for (j in target_marker) {
	    target_cells <- formatted_data[which(formatted_data[,column] == j),     ]
	    target_cell_cords <- target_cells[, c("Cell.ID", "Cell.X.Position","Cell.Y.Position")]
	    target_cell_cords <- remove_rownames(target_cell_cords)
	    target_cell_cords <- target_cell_cords %>% column_to_rownames("Cell.ID")
	    reference_target_result <- frNN(target_cell_cords, eps = radius,
                                         query = reference_cell_cords, sort =FALSE)
	    n_targets <- rapply(reference_target_result$id, length)
        dataframe[,j] <- n_targets
    }
  all.df[[i]] <- dataframe
  }
  return(all.df)
}
