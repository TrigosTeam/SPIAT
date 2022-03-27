#' plot_cell_distances_violin
#'
#' @description Plots distances between cells as a violin plot
#'
#' @param cell_to_cell_dist Data.frame. Output from
#'   \code{\link{calculate_all_distances_between_phenotypes}}.
#' @import ggplot2
#' @examples
#' distances <- calculate_distances_between_cell_types(SPIAT::defined_image,
#' cell_types_of_interest = c("Immune1", "Immune2"), feature_colname="Cell.Type")
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

