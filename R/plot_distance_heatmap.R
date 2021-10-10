#' plot_distance_heatmap
#'
#' @description Takes the output of cell_distances and plot the distances as a heatmap
#'
#' @param phenotype_distances_result Dataframe output from calculate_distances_between_phenotypes
#' @param metric Metric to be plotted. One of "mean", "std.dev" or "median".
#' @import ggplot2
#' @return A plot is returned
#' @examples
#' summary_distances <- calculate_summary_distances_between_cell_types(SPIAT::formatted_image)
#' plot_distance_heatmap(summary_distances)
#' @export

plot_distance_heatmap <- function(phenotype_distances_result, metric = "mean"){
  
    # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
    Reference <- Nearest <- Mean <- Std.Dev <- Median <- NULL

    if(metric == "mean"){
      limit <- range(unlist(phenotype_distances_result$Mean), na.rm=TRUE)
      g <- ggplot(phenotype_distances_result, aes(x = Reference, y = Nearest, fill = Mean)) +
        geom_tile() +
        xlab("Cell of interest (COI)") +
        ylab("Nearest cell to COI") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_fill_viridis_c(limits = limit, direction = -1)
      print(g)

    }else if(metric == "std.dev"){
      limit <- range(unlist(phenotype_distances_result$Std.Dev), na.rm=TRUE)
      g <- ggplot(phenotype_distances_result, aes(x = Reference, y = Nearest, fill = Std.Dev)) + geom_tile() + xlab("Cell of interest (COI)") + ylab("Nearest cell to COI") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_fill_viridis_c(limits = limit, direction = -1)
      print(g)

    }else if(metric == "median"){
      limit <- range(unlist(phenotype_distances_result$Median), na.rm=TRUE)
      g <- ggplot(phenotype_distances_result, aes(x = Reference, y = Nearest, fill = Median)) + geom_tile() + xlab("Cell of interest (COI)") + ylab("Nearest cell to COI") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_fill_viridis_c(limits = limit, direction = -1)
      print(g)
    }
}
