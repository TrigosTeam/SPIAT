#' calculate_summary_distances_between_cell_types
#'
#' @description Returns the mean, median and stardard deviation of the distances
#'   between phenotypes.
#' @param sce_object SingleCellExperiment object in the form of the output of
#'   format_image_to_sce.
#' @param all_combinations Boolean. If TRUE, the distances between all possible
#'   combinations of cell types will be calculated.
#' @param cell_types_of_interest String Vector of cell types to consider if
#'   all_combinations is FALSE.
#' @param feature_colname String. Column of cells to choose the cell type from
#'   (e.g. Cell.Type, Cell.Type2, etc).
#'   
#' @import dplyr
#' @return A data.frame is returned
#' @examples
#' summary_distances <- calculate_summary_distances_between_cell_types(SPIAT::defined_image,
#' feature_colname = "Cell.Type", all_combinations = FALSE, 
#' cell_types_of_interest = c("Tumour","Immune1"))
#' @export

calculate_summary_distances_between_cell_types <- function(sce_object, feature_colname="Cell.Type", 
                                                           all_combinations = FALSE,
                                                           cell_types_of_interest) {
  
  formatted_data <- get_colData(sce_object)
  formatted_data <- formatted_data[,c("Cell.ID","Cell.X.Position", "Cell.Y.Position", feature_colname)]
  formatted_data <- formatted_data[formatted_data[,feature_colname] != "",]
  
  #Get the list of cells under each cell type
  cell_types = list()
  if (all_combinations) {
    # all cell types are cell types of interest
    cell_types_of_interest <- unique(formatted_data[[feature_colname]])
    print("All markers are used in pair-wise distance calculation: ")
    print(unique(formatted_data[[feature_colname]]))
  }else{
    #keep those cells that are selected
    formatted_data <- formatted_data[formatted_data[,feature_colname] %in% cell_types_of_interest, ]
    #CHECK
    if (nrow(formatted_data) == 0){
      stop("No cells belong to the specified marker combinations of interest")
    }
    print("Markers had been selected in pair-wise distance calculation: ")
    print(unique(formatted_data[[feature_colname]]))
  }
  
  # calculate all distances between cell types
  pairwise_dists <- calculate_distances_between_cell_types(sce_object, 
                                         cell_types_of_interest = cell_types_of_interest, 
                                         feature_colname = feature_colname)
  
  # summarise the results
  summarised_dists <- pairwise_dists %>% 
    group_by(Pair) %>%
    summarise(mean(Distance), min(Distance), max(Distance),
              stats::median(Distance), stats::sd(Distance))
    
  summarised_dists <- data.frame(summarised_dists)
  colnames(summarised_dists) <- c("Pair" , "Mean", "Min", "Max", "Median", "Std.Dev")
  
  return(summarised_dists)
}
