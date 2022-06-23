#' calculate_pairwise_distances_between_celltypes
#'
#' @description Returns the pairwise distances between cells of different types.
#'   If none of the cell types are found, it will print an error message and
#'   return a vector of NAs.
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param cell_types_of_interest String Vector containing cell types to be
#'   considered, if NULL, all cell type combinations will be calculated.
#' @param feature_colname String of the name the feature column with the cell
#'   types of interest to be considered.
#' @import dplyr
#' @return A data.frame is returned.
#' @examples
#' calculate_pairwise_distances_between_celltypes(SPIAT::defined_image,
#' cell_types_of_interest = c("Tumour","Immune1"), 
#' feature_colname = "Cell.Type")
#' @export

calculate_pairwise_distances_between_celltypes <- function(
        spe_object, cell_types_of_interest = NULL, feature_colname){
    
    #Reads the image file and deletes cell rows with NA positions
    dat <- get_colData(spe_object)
    
    #Selects all rows in the data file which only contains the cells of interest
    if(!is.null(cell_types_of_interest)){
        unique_cell_types_selected <- cell_types_of_interest
        dat <- dat[dat[, feature_colname] %in% unique_cell_types_selected,]
    }
    
    #CHECK
    if (nrow(dat) == 0) {
        methods::show("There are no cells or no cells of specified cell types")
        cell_to_cell_dist_all <- c(Cell1 = NA, Cell2 = NA, Distance = NA, 
                                   Pair = NA, Type1 = NA, Type2 = NA) 
    }else{
        dat <- dat[,c("Cell.ID", feature_colname, "Cell.X.Position", 
                      "Cell.Y.Position")]
        
        #Creates a list of the number of cell types with all their 
        #corresponding cell ID's
        cell_types <- list()
        for (eachType in unique(dat[, feature_colname])) {
            cell_types[[eachType]] <- as.character(
                dat$Cell.ID[dat[, feature_colname] == eachType])
        }
        
        cell_id_vector <- dat$Cell.ID
        
        #Calculates cell to cell distances
        dist_all <- -apcluster::negDistMat(dat[,c("Cell.X.Position", 
                                                  "Cell.Y.Position")])
        colnames(dist_all) <- cell_id_vector
        rownames(dist_all) <- cell_id_vector
        
        cell_to_cell_dist_all <- vector()
        for(cell_name1 in names(cell_types)){
            for(cell_name2 in names(cell_types)){
                cell_ids1 <- cell_types[[cell_name1]]
                cell_ids2 <- cell_types[[cell_name2]]
                if(length(cell_ids1) >= 3 & length(cell_ids2) >=3){
                    cell_to_cell <- dist_all[cell_id_vector %in% cell_ids1, 
                                             cell_id_vector %in% cell_ids2]
                    #Melts dist_all to produce dataframe of target and nearest 
                    #cell ID's columns and distance column
                    cell_to_cell_dist <- reshape2::melt(cell_to_cell)
                    cell_to_cell_dist$Type1 <- cell_name1
                    cell_to_cell_dist$Type2 <- cell_name2
                    cell_to_cell_dist$Pair <- paste(cell_name1,
                                                    cell_name2,sep="/")
                    
                    if(cell_name1 == cell_name2){
                        cell_to_cell_dist$value[
                            cell_to_cell_dist$value == 0] <- NA
                    }
                    colnames(cell_to_cell_dist)[3] <- "Distance" 
                    cell_to_cell_dist_all <- rbind(cell_to_cell_dist_all, 
                                                   cell_to_cell_dist)
                }
            }
        }
        # remove NAs e.g. for distance of cell against itself
        cell_to_cell_dist_all <- 
            cell_to_cell_dist_all[stats::complete.cases(cell_to_cell_dist_all),]
        colnames(cell_to_cell_dist_all)[c(1,2)] <- c("Cell1", "Cell2")
    }
    
    return(cell_to_cell_dist_all)
}
