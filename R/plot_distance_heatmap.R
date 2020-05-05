#' plot_distance_heatmap
#'
#' @description Takes the output of cell_distances and plot the distances as a heatmap
#'
#' @param phenotype_distances_result Dataframe output from calculate_distances_between_phenotypes
#' @param metric Metric to be plotted. One of "mean", "std.dev" or "median".
#' @import ggplot2
#' @export

plot_distance_heatmap <- function(phenotype_distances_result, metric = "mean"){

    if(metric == "mean"){
      limit <- range(unlist(phenotype_distances_result$Mean), na.rm=TRUE)
      g <- ggplot(phenotype_distances_result, aes(x = Target, y = Nearest, fill = Mean)) +
        geom_tile() +
        xlab("Cell of interest (COI)") +
        ylab("Nearest cell to COI") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), text=element_text(size=14,family="serif")) +
        scale_fill_viridis_c(limits = limit)
      print(g)

    }else if(metric == "std.dev"){
      limit <- range(unlist(phenotype_distances_result$Std.Dev), na.rm=TRUE)
      g <- ggplot(phenotype_distances_result, aes(x = Target, y = Nearest, fill = Std.Dev)) + geom_tile() + xlab("Cell of interest (COI)") + ylab("Nearest cell to COI") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), text=element_text(size=14,family="serif")) +
        scale_fill_viridis_c(limits = limit)
      print(g)

    }else if(metric == "median"){
      limit <- range(unlist(phenotype_distances_result$Median), na.rm=TRUE)
      g <- ggplot(phenotype_distances_result, aes(x = Target, y = Nearest, fill = Median)) + geom_tile() + xlab("Cell of interest (COI)") + ylab("Nearest cell to COI") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), text=element_text(size=14,family="serif")) +
        scale_fill_viridis_c(limits = limit)
      print(g)
    }
}
