#' average_percentage_of_cells_within_radius
#'
#' @description Calculates the average percentage of cells of a target cell type
#'   within a radius from the cells with a reference cell type. The calculation
#'   is done per reference cell, so runtime will depend on the number of
#'   reference cells present. Output is a single value (the mean for the image).
#'
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param reference_celltype String specifying the cell type of reference
#'   cells.
#' @param target_celltype String specifying the cell type for target cells
#' @param radius Integer specifying the radius of search for cells around the
#'   reference cells. Radii of ~100 are recommended. If too small, too few cells
#'   might be present.
#' @param feature_colname String specifying the column with the desired cell
#'   type annotations.
#' @import dplyr
#' @return A numeric vector and a plot are returned
#' @examples
#' average_percentage_of_cells_within_radius(SPIAT::defined_image, "Tumour", 
#' "Immune3", radius = 100, "Cell.Type")
#' @export

average_percentage_of_cells_within_radius <- function(spe_object, 
                                                      reference_celltype, 
                                                      target_celltype, 
                                                      radius = 100, 
                                                      feature_colname){
    # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
    phenotype_names <- output_phenotype <- NULL
    formatted_data <- get_colData(spe_object)
    #Select cells with the reference phenotype
    reference_celltypes <- formatted_data[formatted_data[,feature_colname] 
                                          == reference_celltype,]
    target_celltypes <- formatted_data[formatted_data[,feature_colname] 
                                       == target_celltype,]
    
    #CHECK
    if (nrow(reference_celltypes) == 0 || nrow(target_celltypes) == 0) {
        methods::show("There are no reference cells or no target cells")
        return(NA)
    }else{
        #get the coordinates to find neighbours
        reference_cell_cords <- reference_celltypes[,c("Cell.X.Position", 
                                                       "Cell.Y.Position")]
        
        #frNN output ids
        search_result <- dbscan::frNN(formatted_data[,c("Cell.X.Position", 
                                                        "Cell.Y.Position")],
                                      eps = radius, query = reference_cell_cords, 
                                      sort = FALSE)
        rownums <- unique(unlist(search_result$id))
        
        #CHECK
        if (length(rownums) == 0) {
            methods::show("There are no target cells within the radius")
            return(NA)
        } else {
            output_percentage <- vector()
            cell_IDs <- vector()
            
            for(cell in names(search_result$id)){
                radius_cells <- search_result$id[[cell]]
                radius_cells <- radius_cells[radius_cells != cell]
                radius_cells <- formatted_data[radius_cells,]
                
                percentage <- (nrow(radius_cells[radius_cells$Cell.ID %in% 
                                                     target_celltypes$Cell.ID, ])
                               /nrow(radius_cells))*100
                output_percentage <- c(output_percentage, percentage)
                cell_IDs <- c(cell_IDs, formatted_data[cell,"Cell.ID"])
            }
            names(output_percentage) <- cell_IDs
        }
        return(mean(output_percentage))
    }
}
