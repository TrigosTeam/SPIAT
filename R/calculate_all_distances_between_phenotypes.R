#' calculate_all_distances_between_phenotypes
#'
#' @description Returns the distances between cells of different phenotypes
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param remove_other If TRUE, the plotting will exclude cells with only DAPI marker
#' @param cell_phenotypes_of_interest Vector containing phenotypes to be considered,
#' if NULL, all phenotype combinations will be calculated
#' @import dplyr
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble rownames_to_column
#' @importFrom stats complete.cases
#' @importFrom apcluster negDistMat
#' @importFrom reshape2 melt
#' @export

calculate_all_distances_between_phenotypes <- function(sce_object, remove_other = TRUE, cell_phenotypes_of_interest = NULL){

    #Reads the image file and deletes cell rows with NA positions
    dat <- data.frame(colData(sce_object))
    dat <- dat[complete.cases(dat),]
    dat<- dat %>% rownames_to_column("Cell.ID") #convert rowname to column
    
    #Selects all rows in the data file which only contains the cells of interest
    if(!is.null(cell_phenotypes_of_interest)){
        unique_cell_phenotypes_selected <- as.vector(unique(unlist(cell_phenotypes_of_interest)))
        dat <- dat[dat$Phenotype %in% unique_cell_phenotypes_selected,]
    }
    #CHECK
    if (nrow(dat) == 0) {
      stop("There are no cells or no cells of specified phenotypes")
    }

    dat <- dat[,c("Cell.ID","Phenotype", "Cell.X.Position", "Cell.Y.Position")]
    dat <- dat[dat$Phenotype != "",]
    dat$Phenotype <- as.character(dat$Phenotype)
    if (isTRUE(remove_other)){
        dat <- dat[dat$Phenotype != "OTHER",]
    }

    #Creates a list of the number of cell types with all their corresponding cell ID's
    cell_phenotypes = list()
    for (eachType in unique(dat$Phenotype)) {
        cell_phenotypes[[eachType]] = as.character(dat$Cell.ID[dat$Phenotype == eachType])
    }

    cell_id_vector <- dat$Cell.ID

    #Calculates cell to cell distances
    dist_all <- - negDistMat(dat[,c("Cell.X.Position", "Cell.Y.Position")])
    colnames(dist_all) <- cell_id_vector
    rownames(dist_all) <- cell_id_vector

    cell_to_cell_dist_all <- vector()
    for(cell_name1 in names(cell_phenotypes)){
      for(cell_name2 in names(cell_phenotypes)){
          cell_ids1 <- cell_phenotypes[[cell_name1]]
          cell_ids2 <- cell_phenotypes[[cell_name2]]

          if(length(cell_ids1) >= 3 & length(cell_ids2) >=3){
            cell_to_cell <- dist_all[cell_id_vector %in% cell_ids1, cell_id_vector %in% cell_ids2]
            
            #Melts dist_all to produce dataframe of target and nearest cell ID's columns and distance column
            cell_to_cell_dist <- melt(cell_to_cell)
            cell_to_cell_dist$Pair <- paste(cell_name1,cell_name2,sep="_")
            
            if(cell_name1 == cell_name2){
              cell_to_cell_dist$value[cell_to_cell_dist$value == 0] <- NA
            }
            colnames(cell_to_cell_dist)[3] <- "Distance" 
            cell_to_cell_dist_all <- rbind(cell_to_cell_dist_all, cell_to_cell_dist)
          }
      }
    }
  
    # remove NAs e.g. for distance of cell against itself
    cell_to_cell_dist_all <- cell_to_cell_dist_all[complete.cases(cell_to_cell_dist_all),]
    
  return(cell_to_cell_dist_all)
}

