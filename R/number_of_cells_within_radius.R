#' Number of cells within a radius
#'
#' @description Calculates the number of cells of a target cell type
#'   within a pre-defined radius around cells of a reference cell type.
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param reference_celltype String. Cell type to be used for reference cells.
#' @param target_celltype String. Cell type to be used for target cells.
#' @param radius Numeric. Radius around the reference cells.
#' @param feature_colname String specifying the column with the desired cell
#'   type annotations.
#' @import dplyr
#' @return A list of dataframes with the number of target cells of each of the
#'   reference cells
#' @export
#' @examples
#' n_in_radius <- number_of_cells_within_radius(SPIAT::defined_image,
#' reference_celltype = "Tumour", target_celltype="Immune1", radius = 50,
#' feature_colname = "Cell.Type")

number_of_cells_within_radius <- function(spe_object, reference_celltype, 
                                           target_celltype, radius = 20, 
                                           feature_colname) 
{
    formatted_data <- get_colData(spe_object)
    
    all.df <- list()
    for (i in reference_celltype) {
        reference_cells <- formatted_data[which(formatted_data[,feature_colname] == i), ]
        reference_cell_cords <- reference_cells[, c("Cell.ID", "Cell.X.Position","Cell.Y.Position")]
        dataframe <- tibble::remove_rownames(reference_cell_cords)
        dataframe <- dataframe %>% tibble::column_to_rownames("Cell.ID")
        reference_cell_cords <- reference_cells[, c( "Cell.X.Position","Cell.Y.Position")]
        for (j in target_celltype) {
            target_cells <- formatted_data[which(formatted_data[,feature_colname] == j),     ]
            target_cell_cords <- target_cells[, c("Cell.ID", "Cell.X.Position","Cell.Y.Position")]
            target_cell_cords <- tibble::remove_rownames(target_cell_cords)
            target_cell_cords <- target_cell_cords %>% tibble::column_to_rownames("Cell.ID")
            reference_target_result <- dbscan::frNN(target_cell_cords, eps = radius,
                                                    query = reference_cell_cords, sort =FALSE)
            n_targets <- rapply(reference_target_result$id, length)
            dataframe[,j] <- n_targets
        }
        all.df[[i]] <- dataframe
    }
    return(all.df)
}
