#' average_minimum_distance
#'
#' @description Calculates the average minimum distance of all cells to their
#'   nearest cells in the input image.
#'
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @import dplyr
#' @return A single number is returned
#' @examples
#' average_minimum_distance(SPIAT::simulated_image)
#' @export

average_minimum_distance <- function(spe_object) {

    formatted_data <- get_colData(spe_object)

    #extract the cell coordinates
    all_cell_cords <- formatted_data[,c("Cell.X.Position", "Cell.Y.Position")]
    
    #CHECK
    if (nrow(all_cell_cords) == 0) {
        stop("No cells found in average minimum distance calculation")
    }
    
    #calculate the closest 2 neighbours, 1st being itself
    all_closest <- RANN::nn2(data = all_cell_cords, k = 2)

    #grab the distances and find the average
    all_closest_dist <- all_closest$nn.dists[,2]
    average_min_distance <- mean(all_closest_dist)

    return(average_min_distance)
}
