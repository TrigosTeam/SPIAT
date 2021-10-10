#' calculate_minimum_distances
#'
#' @description Returns the distance of the closest cell of a specific type from each reference cell
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param cell_types_of_interest Vector of marker combinations to consider
#' is FALSE
#' @param column Column of cells to choose the phenotype from (e.g. Cell.Type, Cell.Type2, etc)
#' @importFrom apcluster negDistMat
#' @importFrom RANN nn2
#' @importFrom gtools permutations
#' @import dplyr
#' @importFrom stats median sd
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble rownames_to_column
#' @return A data.frame is returned
#' @examples
#' @export

calculate_minimum_distances <- function(sce_object, column="Phenotype",
                                                           cell_types_of_interest = NULL) {
  
  formatted_data <- data.frame(colData(sce_object))
  formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column
  
  formatted_data <- formatted_data[,c("Cell.ID","Cell.X.Position", "Cell.Y.Position", column)]
  formatted_data <- formatted_data[formatted_data[,column] != "",]
  
  #Add a new column Cell_type which duplicates the phenotype
  formatted_data$Cell_type <- as.character(formatted_data[,column])
  
  #Get the list of cells under each cell type
  cell_types = list()
  #assign cells to each specified cell type
  for (eachType in cell_types_of_interest) {
      cell_types[[eachType]] = as.character(formatted_data$Cell.ID[formatted_data[,column] == eachType])
    }
  #keep those cells that are selected
  formatted_data <- formatted_data[formatted_data[,column] %in% cell_types_of_interest, ]
  #CHECK
  if (nrow(formatted_data) == 0){
      stop("No cells belong to the specified marker combinations of interest")
    }
  print("Markers had been selected in pair-wise distance calculation: ")
  print(unique(formatted_data$Cell_type))
  
  #permutation for the different cell type combinations
  permu = permutations(length(unique(formatted_data$Cell_type)), 2, repeats.allowed = TRUE)
  unique_cells <- unique(formatted_data$Cell_type) #unique cell types
  result = vector()
  
  for (i in seq_len(nrow(permu))) {
    eachPermu = permu[i, ]
    name1 = unique_cells[eachPermu[1]]
    name2 = unique_cells[eachPermu[2]]
    
    cell_type1 <- as.character(cell_types[[name1]])
    cell_type2 <- as.character(cell_types[[name2]])
    
    #No need to calculate min distance when they're the same cell type
    if (name1 == name2) {
       local_dist_mins <- data.frame(RefCell=NA,
                                    RefType = name1,
                                    NearestCell = NA,
                                    NearestType = name2,
                                    Dist = NA)
    } else {
      #vector to store all mins
      local_dist_min <- vector()
      #grab coordinates of all other cells in cell_type2 to calculate distances with
      all_celltype2_cord <- formatted_data[formatted_data[,column] == name2, c("Cell.X.Position", "Cell.Y.Position")]
      all_celltype1_cord <- formatted_data[formatted_data[,column] == name1, c("Cell.X.Position", "Cell.Y.Position")]
      
      #find all of closest points
      all_closest <- nn2(data = all_celltype2_cord, query = all_celltype1_cord, k = 1)
      
      all_celltype2_cord2 <- formatted_data[formatted_data[,column] == name2, c("Cell.ID", "Cell.X.Position", "Cell.Y.Position")]
      
      local_dist_mins <- all_closest$nn.dists
      local_dist_mins <- as.vector(local_dist_mins)
      local_dist_mins <- data.frame(RefCell=formatted_data[formatted_data[,column] == name1, c("Cell.ID")],
                                    RefType = name1,
                                    NearestCell = all_celltype2_cord2$Cell.ID[as.vector(all_closest$nn.idx)],
                                    NearestType = name2,
                                    Dist = local_dist_mins)
    }
    result <- rbind(result, local_dist_mins)
  }
  
  # remove NAs e.g. for distance of cell against itself
  result <- result[complete.cases(result),]
  
  return(result)
}