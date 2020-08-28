#' calculate_cell_proportions
#'
#' @description Calculate the number and proportion of each cell phenotype in image
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param markers_of_interest Vector specifying markers to calculate proportions for (optional)
#' @importFrom SummarizedExperiment colData
#' @importFrom stats complete.cases
#' @export

calculate_cell_proportions <- function(sce_object, markers_of_interest = NULL, markers_to_exclude = NULL){

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
    
    if (!is.null(markers_to_exclude)) {
        cell_proportions$Proportion[cell_proportions$Cell_type %in% markers_to_exclude] <- NA
    }
    cell_proportions$Percentage <- cell_proportions$Proportion*100
    cell_proportions$Proportion_name <- paste(cell_proportions$Cell_type, "Total", sep="/")
    
    if (!is.null(markers_of_interest)) {
        
        datalist = list()
        
        for (marker in markers_of_interest) {
            
            marker_cells <- cell_proportions[grepl(marker, cell_proportions$Cell_type), ]
            marker_cells_total <- sum(marker_cells$Number_of_cells)
            
            marker_cells$Reference <- rep(marker, nrow(marker_cells))
            marker_cells$Proportion <- marker_cells$Number_of_cells/marker_cells_total
            marker_cells$Percentage <- marker_cells$Proportion*100
            marker_cells$Proportion_name <- paste(marker_cells$Cell_type, marker, sep="/")
            
            datalist[[marker]] <- marker_cells
            
        }
        
        cell_proportions <- do.call(rbind, datalist)
    }
    
    cell_proportions <- cell_proportions[rev(order(cell_proportions$Proportion)), ]
    rownames(cell_proportions) <- NULL
    
    return(cell_proportions)
}
