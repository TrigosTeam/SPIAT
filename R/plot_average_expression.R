#' plot_average_expression
#'
#' @description Takes in a vector or radii and calculates the average expression
#' of a target marker using average_expression function. It plots the expression
#' level as a line graph
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param reference_marker String specifying the reference marker
#' @param target_marker String specifying the marker to calculate its average expression
#' @param radii Vector of integers specifying the search radius around reference cells
#' @export

plot_average_expression <- function(sce_object, reference_marker, target_marker, radii) {

    average_expression_result <- vector()

    for (radius in radii) {
        result <- average_marker_expression_within_radius(sce_object, reference_marker, target_marker, radius = radius)
        #check
        if(!is.numeric(result)) {
            return(print(paste("Cannot calculate average expression for radius = ", radius, sep = "")))
        }
        average_expression_result <- c(average_expression_result, result)
    }

    average_expression_df <- data.frame(radii, average_expression_result)

    p <- ggplot(data=average_expression_df, aes(x=radii, y=average_expression_result))

    title <- paste("Average expression of ", target_marker, " from ", reference_marker, " for various radii", sep="")
    p <- p + geom_line() + geom_point() + ggtitle(title) + xlab("Radius") + ylab("Average Expression")
    p
}
