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
    #assign cells to each cell type based on the markers expressed
    for (eachType in unique(formatted_data[[feature_colname]])) {
      cell_types[[eachType]] = as.character(formatted_data$Cell.ID[formatted_data[,feature_colname] == eachType])
    }
    print("All markers are used in pair-wise distance calculation: ")
    print(unique(formatted_data[[feature_colname]]))
  }else{
    #assign cells to each specified cell type
    for (eachType in cell_types_of_interest) {
      cell_types[[eachType]] = as.character(formatted_data$Cell.ID[formatted_data[,feature_colname] == eachType])
    }
    #keep those cells that are selected
    formatted_data <- formatted_data[formatted_data[,feature_colname] %in% cell_types_of_interest, ]
    #CHECK
    if (nrow(formatted_data) == 0){
      stop("No cells belong to the specified marker combinations of interest")
    }
    print("Markers had been selected in pair-wise distance calculation: ")
    print(unique(formatted_data[[feature_colname]]))
  }
  
  #permutation for the different cell type combinations
  comb = gtools::permutations(length(unique(formatted_data[[feature_colname]])), 2, repeats.allowed = T)
  unique_cells <- unique(formatted_data[[feature_colname]]) #unique cell types
  result = vector()
  
  for (i in seq_len(nrow(comb))) {
    eachComb = comb[i, ]
    name1 = unique_cells[eachComb[1]]
    name2 = unique_cells[eachComb[2]]
    
    cell_type1 <- as.character(cell_types[[name1]])
    cell_type2 <- as.character(cell_types[[name2]])
    
    #No need to calculate min distance when they're the same cell type
    if (name1 == name2) {
      local_result = data.frame(Reference = name1, Nearest = name2, Mean = NA, Std.Dev = NA, Median = NA)
    } else {
      #vector to store all mins
      local_dist_min <- vector()
      #grab coordinates of all other cells in cell_type2 to calculate distances with
      all_celltype2_cord <- formatted_data[formatted_data[,feature_colname] == name2, c("Cell.X.Position", "Cell.Y.Position")]
      all_celltype1_cord <- formatted_data[formatted_data[,feature_colname] == name1, c("Cell.X.Position", "Cell.Y.Position")]
      
      #find all of closest points
      all_closest <- RANN::nn2(data = all_celltype2_cord, query = all_celltype1_cord, k = 1)
      
      local_dist_mins <- all_closest$nn.dists
      
      local_result <- data.frame(Reference = name1, Nearest = name2, 
                                 Mean = mean(local_dist_mins), 
                                 Std.Dev = stats::sd(local_dist_mins), 
                                 Median = stats::median(local_dist_mins))
    }
    result <- rbind(result, local_result)
  }
  
  # remove NAs e.g. for distance of cell against itself
  result <- result[complete.cases(result),]
  
  return(result)
}
