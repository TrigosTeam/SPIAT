#' percentage_of_cells_within_radius
#'
#' @description Calculates the percentage of cells of a target phenotype within a radius
#' from the cells with a reference phenotype.
#' The calculation is done per reference cell, so runtime will depend on the number of reference
#' cells present.
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param reference_phenotypes String specifying the phenotypes of reference cells
#' @param target_phenotypes String specifying the phenotypes for target cells
#' @param radius Integer specifying the radius of search for cells around the reference cells.
#' Radii of ~100 are recommended. If too small, too few cells might be present
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom dbscan frNN
#' @import stats
#' @import SingleCellExperiment
#' @export

# %>% operator is in package 'magrittr' but imported by dplyr
# colData() is in package 'SummarizedExperiment' but imported by SingleCellExperiment

percentage_of_cells_within_radius <- function(sce_object, reference_phenotypes, target_phenotypes, radius = 100){

  formatted_data <- data.frame(colData(sce_object))
  formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

  #Select cells with the reference phenotype
  reference_phenotypes <- formatted_data[formatted_data$Phenotype == reference_phenotypes,]
  target_phenotypes <- formatted_data[formatted_data$Phenotype == target_phenotypes,]
  
  #CHECK
  if (nrow(reference_phenotypes) == 0 || nrow(target_phenotypes) == 0) {
    stop("There are no reference cells or target cells, calculation aborted")
  }

  #get the coordinates to find neighbours
  reference_cell_cords <- reference_phenotypes[,c("Cell.X.Position", "Cell.Y.Position")]

  #frNN output ids
  search_result <- frNN(formatted_data[,c("Cell.X.Position", "Cell.Y.Position")],
                        eps = radius, query = reference_cell_cords, sort = FALSE)

  rownums <- unique(unlist(search_result$id))

  #CHECK
  if (length(rownums) == 0) {
    stop("There are no target cells within the radius")
  } else {
    output_percentage <- vector()
    cell_IDs <- vector()

    for(cell in names(search_result$id)){
      radius_cells <- search_result$id[[cell]]
      radius_cells <- radius_cells[radius_cells != cell]
      radius_cells <- formatted_data[radius_cells,]
      percentage <- (sum(radius_cells$Phenotype == target_phenotypes)/nrow(radius_cells))*100
      output_percentage <- c(output_percentage, percentage)
      cell_IDs <- c(cell_IDs, formatted_data[cell,"Cell.ID"])
    }
    names(output_percentage) <- cell_IDs
  }

  return(output_percentage)
}
