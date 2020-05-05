#' average_minimum_distance
#'
#' @description Calculates the average minimum distance of all cells in the sce_object
#'
#' @param sce_object Singlecellexperiment object in the form of the output of format_image_to_sce
#' @importFrom RANN nn2
#' @import dplyr
#' @import SingleCellExperiment
#' @importFrom tibble rownames_to_column
#' @export

# %>% operator is in package 'magrittr' but imported by dplyr
# colData() is in package 'SummarizedExperiment' but imported by SingleCellExperiment

average_minimum_distance <- function(sce_object) {

    formatted_data <- data.frame(colData(sce_object))
    #formatted_data <- formatted_data[complete.cases(formatted_data),]
    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

    #extract the cell coordinates
    all_cell_cords <- formatted_data[,c("Cell.X.Position", "Cell.Y.Position")]
    
    #CHECK
    if (nrow(all_cell_cords) == 0) {
        stop("No cells found in average minimum distance calculation")
    }
    
    #calculate the closest 2 neighbours, 1st being itself
    all_closest <- nn2(data = all_cell_cords, k = 2)

    #grab the distances and find the average
    all_closest_dist <- all_closest$nn.dists[,2]
    average_min_distance <- mean(all_closest_dist)

    return(average_min_distance)
}
