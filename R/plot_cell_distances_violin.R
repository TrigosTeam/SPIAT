#' plot_cell_distances_violin
#'
#' @description Plots distances between cells as a violin plot
#'
#' @param cell_to_cell_dist Output from calculate_all_distances_between_phenotypes
#' @import ggplot2
#' @examples
#' distances <- calculate_all_distances_between_cell_types(SPIAT::formatted_image, cell_types_of_interest = c("CD3,CD4", "CD3,CD8"), column="Phenotype")
#' plot_cell_distances_violin(distances)
#' @return A plot is returned
#' @export


plot_cell_distances_violin <- function(cell_to_cell_dist){
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Pair <- Distance <- NULL
  
  ggplot(cell_to_cell_dist, aes(x = Pair, y = Distance)) + geom_violin() +
    facet_wrap(~Pair, scales="free_x") +
    theme_bw()
  
}

