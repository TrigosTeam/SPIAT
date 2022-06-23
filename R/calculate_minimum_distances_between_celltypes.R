#' calculate_minimum_distances_between_celltypes
#'
#' @description Returns the distance of the closest cell of a specific type from
#'   each reference cell.
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param cell_types_of_interest String Vector of marker combinations to
#'   consider is FALSE.
#' @param feature_colname String of the feature column of cells to choose the
#'   cell types from (e.g. Cell.Type, Cell.Type2, etc).
#' @import dplyr
#' @return A data.frame is returned
#' @examples
#' min_dists <- calculate_minimum_distances_between_celltypes(
#' SPIAT::defined_image, feature_colname = "Cell.Type", 
#' cell_types_of_interest = c("Tumour","Immune1"))
#' @export

calculate_minimum_distances_between_celltypes <- function(spe_object, 
                                                          feature_colname,
                                                          cell_types_of_interest 
                                                          = NULL) {
    Pair <- Distance <- NULL
    formatted_data <- get_colData(spe_object)
    
    formatted_data <- formatted_data[,c("Cell.ID","Cell.X.Position", 
                                        "Cell.Y.Position", feature_colname)]
    formatted_data <- formatted_data[formatted_data[,feature_colname] != "",]
    
    #Get the list of cells under each cell type
    cell_types <- list()
    #assign cells to each specified cell type
    for (eachType in cell_types_of_interest) {
        cell_types[[eachType]] <- as.character(formatted_data$Cell.ID[
            formatted_data[,feature_colname] == eachType])
    }
    #keep those cells that are selected
    formatted_data <- formatted_data[formatted_data[,feature_colname] %in% 
                                         cell_types_of_interest, ]
    #CHECK
    if (nrow(formatted_data) == 0){
        stop("No cells belong to the specified marker combinations of interest")
    }
    methods::show("Markers had been selected in minimum distance calculation: ")
    methods::show(unique(formatted_data[[feature_colname]]))
    
    #different cell type combinations
    permu <- gtools::permutations(length(unique(
        formatted_data[[feature_colname]])), 2, repeats.allowed = TRUE)
    unique_cells <- unique(formatted_data[[feature_colname]]) #unique cell types
    result <- vector()
    
    for (i in seq_len(nrow(permu))) {
        eachPermu <- permu[i, ]
        name1 <- unique_cells[eachPermu[1]]
        name2 <- unique_cells[eachPermu[2]]
        
        cell_type1 <- as.character(cell_types[[name1]])
        cell_type2 <- as.character(cell_types[[name2]])
        
        #No need to calculate min distance when they're the same cell type
        if (name1 == name2) {
            local_dist_mins <- data.frame(RefCell=NA,
                                          RefType = name1,
                                          NearestCell = NA,
                                          NearestType = name2,
                                          Distance = NA)
        } else {
            #vector to store all mins
            local_dist_min <- vector()
            #grab coordinates of all other cells in cell_type2 to calculate å
            all_celltype2_cord <- 
                formatted_data[formatted_data[,feature_colname] == name2, 
                               c("Cell.X.Position", "Cell.Y.Position")]
            all_celltype1_cord <- 
                formatted_data[formatted_data[,feature_colname] == name1, 
                               c("Cell.X.Position", "Cell.Y.Position")]
            
            #find all of closest points
            all_closest <- RANN::nn2(data = all_celltype2_cord, 
                                     query = all_celltype1_cord, k = 1)
            
            all_celltype2_cord2 <- 
                formatted_data[formatted_data[,feature_colname] == name2, 
                               c("Cell.ID", "Cell.X.Position", 
                                 "Cell.Y.Position")]
            
            local_dist_mins <- all_closest$nn.dists
            local_dist_mins <- as.vector(local_dist_mins)
            local_dist_mins <- data.frame(
                RefCell=formatted_data[formatted_data[,feature_colname] == 
                                           name1, c("Cell.ID")],
                RefType = name1,
                NearestCell = 
                    all_celltype2_cord2$Cell.ID[as.vector(all_closest$nn.idx)],
                NearestType = name2,
                Distance = local_dist_mins)
        }
        result <- rbind(result, local_dist_mins)
    }
    
    # remove NAs e.g. for distance of cell against itself
    result <- result[stats::complete.cases(result),]
    result$Pair <- paste(result$RefType, result$NearestType,sep = "/")
    
    return(result)
}
