#' plot_cell_marker_levels
#'
#' @description Produces a scatter plot of the level of a marker in each cell.
#'   The level of the marker in all cells is shown, at x-y positions, no matter
#'   if cells are phenotyped as being positive or negative for the particular
#'   marker.
#'
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param marker String. Marker to plot.
#' @import dplyr
#' @import ggplot2
#' @return A plot is returned
#' @examples
#' plot_cell_marker_levels(SPIAT::simulated_image, "Immune_marker1")
#' @export

plot_cell_marker_levels <- function(spe_object, marker) {

    Cell.X.Position <- Cell.Y.Position <- NULL
    
    intensity_matrix <- SummarizedExperiment::assay(spe_object)
    markers <- rownames(intensity_matrix)
    
    #CHECK
    if (is.element(marker, markers) == FALSE) {
        stop("The marker specified is not in the data")
    }
    
    formatted_data <- bind_info(spe_object)
    
    #selecting cells that do not contain the marker
    #for one entry that is not marker
    rows <- formatted_data[formatted_data$Phenotype != marker, ] 
    #for multiple entries that does not contain marker
    rows <- rows[!grepl(marker, rows$Phenotype), ] 
    
    #for those cell without the marker, set marker intensity to 0
    #and merge the formatted_data
    rows[, marker] <- 0
    formatted_data[match(rows$Cell.ID,formatted_data$Cell.ID),]<-rows
    
    #selecting the cells that have intensity for a specific marker
    column <- which(colnames(formatted_data) == marker)
    rows_non_zero <- which(formatted_data[,column] != 0)
    intensity_by_marker <- formatted_data[rows_non_zero,]
    
    if (nrow(intensity_by_marker) == 0) {
        methods::show(paste("There are no true intensity for: ", marker, sep=""))
    }
    
    #log the intensity to improve contrast
    intensity_by_marker[,marker] <- log10(intensity_by_marker[,marker])
    
    ggplot(intensity_by_marker, aes(x = Cell.X.Position, y = Cell.Y.Position,
                                    colour = eval(parse(text = marker)))) +
        geom_point(aes(colour=eval(parse(text = marker))),size = 0.1) +
        ggtitle(marker) +
        guides(alpha = "none") + scale_colour_viridis_c(direction = -1) +
        labs(colour = paste("log10","(", as.character(marker),
                            " Intensity", ")",sep="")) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white"),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
}
