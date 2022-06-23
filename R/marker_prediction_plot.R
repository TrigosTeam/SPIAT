#' marker_prediction_plot
#'
#' @description Takes in the returned dataframe from marker_threshold_plot and
#' generates a .pdf file containing scatter plots of actual intensity and
#' predicted intensity for every marker.
#'
#' @param predicted_data Output from \code{\link{predict_phenotypes}}.
#' @param marker String. Marker to plot
#' @import ggplot2
#' @return A plot is returned
#' @examples
#' predicted_result <- predict_phenotypes(spe_object = simulated_image, thresholds = NULL,
#' tumour_marker = "Tumour_marker",baseline_markers = c("Immune_marker1", "Immune_marker2", 
#' "Immune_marker3", "Immune_marker4"), reference_phenotypes = TRUE)
#' marker_prediction_plot(predicted_result, marker = "Tumour_marker")
#' @export

marker_prediction_plot <- function(predicted_data, marker) {
    
    # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
    Cell.X.Position <- Cell.Y.Position <- NULL
    
    #get the markers
    actual_phenotype_colnames <- predicted_data[grepl("_actual_phenotype", colnames(predicted_data))]
    markers <- gsub("_actual_phenotype", "", actual_phenotype_colnames)
    
    #extract the actual and predicted intensity status along with X and Y cord
    actual_colname <- paste(marker, "_actual_phenotype", sep="")
    predicted_colname <- paste(marker, "_predicted_phenotype", sep="")
    
    actual <- predicted_data[,c("Cell.X.Position", "Cell.Y.Position", actual_colname)]
    pred <- predicted_data[,c("Cell.X.Position", "Cell.Y.Position", predicted_colname)]
    
    #plot actual intensity status
    title_actual <- paste("Actual intensity status of ", marker, sep="")
    p_actual <- ggplot(actual, aes(x = Cell.X.Position, y = Cell.Y.Position, color = as.character(actual[,3]))) +
        geom_point(size = 0.1) + scale_color_manual(values=c('grey','red')) +
        guides(alpha = "none") + labs(colour = "Intensity status") + ggtitle(title_actual) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white"),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    
    #plot predicted intensity status
    title_pred <- paste("Predicted intensity status of ", marker, sep="")
    p_pred <- ggplot(pred, aes(x = Cell.X.Position, y = Cell.Y.Position, color = as.character(pred[,3]))) +
        geom_point(size = 0.1) + scale_color_manual(values=c('grey','red')) +
        guides(alpha = "none") + labs(colour = "Intensity status") + ggtitle(title_pred) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white"),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    gridExtra::grid.arrange(p_actual, p_pred, nrow=1)
}
