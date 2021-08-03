#' plot_cell_marker_levels
#'
#' @description Produces a scatter plot of the level of a marker in each cell.
#' The level of the marker in all cells is shown, whether phenotyped as being positive or negative for the particular marker.
#'
#' @param sce_object Singlecellexperiment object in the form of the output of format_image_to_sce
#' @param marker Marker to plot
#' @import dplyr
#' @import ggplot2
#' @importFrom SummarizedExperiment colData assay
#' @importFrom tibble rownames_to_column
#' @return A plot is returned
#' @examples
#' plot_cell_marker_levels(SPIAT::formatted_image, "CD3")
#' @export

plot_cell_marker_levels <- function(sce_object, marker) {
  
    # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
    Cell.X.Position <- Cell.Y.Position <- NULL

    formatted_data <- data.frame(colData(sce_object))

    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

    intensity_matrix <- assay(sce_object)

    markers <- rownames(intensity_matrix)
    
    #CHECK
    if (is.element(marker, markers) == FALSE) {
      stop("The marker specified is not in the data")
    }
    
    cell_ids <- colnames(intensity_matrix)

    rownames(intensity_matrix) <- NULL
    colnames(intensity_matrix) <- NULL
    intensity_matrix_t <- t(intensity_matrix)
    intensity_df <- data.frame(intensity_matrix_t)
    colnames(intensity_df) <- markers

    formatted_data <- cbind(formatted_data, intensity_df)
    formatted_data <- formatted_data[complete.cases(formatted_data),]

    #selecting cells that do not contain the marker
    rows <- formatted_data[formatted_data$Phenotype != marker, ] #for one entry that is not marker
    rows <- rows[!grepl(marker, rows$Phenotype), ] #for multiple entries that does not contain marker
      
    #for those cell without the marker, set marker intensity to 0
    #and merge the formatted_data
    rows[, marker] <- 0
    formatted_data[match(rows$Cell.ID,formatted_data$Cell.ID),]<-rows
        
    #selecting the cells that have intensity for a specific marker
    column <- which(colnames(formatted_data) == marker)
    rows_non_zero <- which(formatted_data[,column] != 0)
    intensity_by_marker <- formatted_data[rows_non_zero,]
        
    if (nrow(intensity_by_marker) == 0) {
      print(paste("There are no true intensity for: ", marker, sep=""))
    }
        
    #log the intensity to improve contrast
    intensity_by_marker[,marker] <- log10(intensity_by_marker[,marker])
    #print(intensity_by_marker)
        
    ggplot(intensity_by_marker, aes(x = Cell.X.Position, y = Cell.Y.Position, colour = eval(parse(text = marker)))) +
      geom_point(aes(colour=eval(parse(text = marker))),size = 0.1) +
      ggtitle(marker) +
      guides(alpha = "none") + scale_colour_viridis_c(direction = -1) +
      labs(colour = paste("log10","(", as.character(marker)," Intensity", ")", sep="")) +
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

