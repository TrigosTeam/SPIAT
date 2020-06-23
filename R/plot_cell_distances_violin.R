#' plot_cell_distances_violin
#'
#' @description Plots distances between cells as a violin plot
#'
#' @param cell_to_cell_dist Output from calculate_all_distances_between_phenotypes
#' @import ggplot2
#' @export


plot_cell_distances_violin <- function(cell_to_cell_dist){
  for(pair in unique(cell_to_cell_dist$Pair)){
    temp <- cell_to_cell_dist[cell_to_cell_dist$Pair == pair,]
    violin_plot <- ggplot(temp, aes(x = Pair, y = Distance)) + geom_violin() +
      ggtitle(paste("Distance between", pair))
    print(violin_plot)
  }
}

