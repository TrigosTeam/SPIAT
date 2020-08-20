#' calculate_cell_proportions
#'
#' @description Calculate the number and proportion of each cell phenotype in image
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @importFrom SummarizedExperiment colData
#' @importFrom stats complete.cases
#' @export

calculate_cell_proportions <- function(sce_object){

    #Reads the image file and deletes cell rows with NA positions
    cell_loc <- data.frame(colData(sce_object))
    cell_loc <- cell_loc[complete.cases(cell_loc),]
    
    #CHECK
    if (nrow(cell_loc) == 0) {
        stop("No cells found for calculating cell proportions")
    }

    #Creates frequency/bar plot of all cell types in the entire image
    cell_proportions <- as.data.frame(table(cell_loc$Phenotype))
    names(cell_proportions)[1] <- 'Cell_type'
    names(cell_proportions)[2] <- 'Number_of_cells'
    total_cells <- sum(cell_proportions$Number_of_cells[cell_proportions$Cell_type != "OTHER"])
    total_cell_proportions <- (cell_proportions$Number_of_cells/total_cells)
    cell_proportions$Proportion <- total_cell_proportions
    cell_proportions$Proportion[cell_proportions$Cell_type == "OTHER"] <- NA
    cell_proportions$Percentage <- cell_proportions$Proportion*100

    return(cell_proportions)
}
