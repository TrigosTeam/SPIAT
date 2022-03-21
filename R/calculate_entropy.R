#' calculate_entropy
#'
#' @description Returns a dataframe of entropies for all reference cells or an
#'   entropy number for the whole image if a radius is not supplied.
#' @param sce_object SingleCellExperiment object in the form of the output of
#'   format_image_to_sce.
#' @param cell_types_of_interest String Vector. Cell types of interest, the
#'   first cell type is considered as reference cell type.
#' @param feature_colname String specifying the column the cell types are from.
#' @param radius (OPTIONAL) Numeric. The maximum radius around a reference cell
#'   for another cell to be considered an interaction.
#' @importFrom SummarizedExperiment colData
#' @return A dataframe or a number depending on the argument radius
#' @export
#' @examples 
#' calculate_entropy(SPIAT::defined_image, cell_types_of_interest = c("Immune1","Immune2"), 
#' feature_colname = "Cell.Type")

calculate_entropy <- function(sce_object, cell_types_of_interest, 
                              feature_colname = "Phenotype", radius = NULL){
  
  if((cell_types_of_interest[1] %in% colData(sce_object)[,feature_colname]) && 
     any(cell_types_of_interest[-1] %in% colData(sce_object)[,feature_colname])){
    if (!is.null(radius)){
      reference_marker <- cell_types_of_interest[1]
      target_marker<-cell_types_of_interest
      n_cells <- number_of_cells_within_radius(sce_object, reference_celltype = reference_marker,
                                               target_celltype = target_marker,radius = radius, feature_colname = feature_colname)
      n_cells.df <- n_cells[[reference_marker]]
      n_cells.df[,"total"] <- 0
      
      for (i in target_marker){
        n_cells.df[,paste(i,"_log2",sep = "")] <- log2(n_cells.df[,i])
        n_cells.df[,"total"] <- n_cells.df[,"total"]+n_cells.df[,i]
      }
      n_cells.df[,"total_log2"] <- log2(n_cells.df[,"total"])
      n_cells.df[n_cells.df == -Inf] <- 0
      
      for (i in target_marker){
        n_cells.df[which(n_cells.df[,"total"] != 0), paste(i,"ratio",sep = "")] <- n_cells.df[,i]/n_cells.df[,"total"]
        n_cells.df[which(n_cells.df[,"total"] == 0), paste(i,"ratio",sep = "")] <- 0
        n_cells.df[,paste(i,"_entropy",sep = "")] <- n_cells.df[, paste(i,"ratio",sep = "")] *(n_cells.df[,paste(i,"_log2",sep = "")] - n_cells.df[,"total_log2"])
      }
      n_cells.df[,"entropy"] <- -(rowSums(n_cells.df[,grepl("_entropy",colnames(n_cells.df))]))  
      
      return(n_cells.df)
    } else{
      entropy_all <- 0
      data <- data.frame(colData(sce_object))
      data <- data[which(data[,feature_colname] %in% cell_types_of_interest),]
      n_all <- dim(data)[1]
      
      if (n_all == 0){
        return(0.0)
      } else{
        for (type in cell_types_of_interest){
          type_data <- data[which(data[,feature_colname]==type),]
          n_type <- dim(type_data)[1]
          if (n_type != 0){
            p_type <- n_type/n_all
            entropy_type <- -(p_type)*log2(p_type)
            entropy_all <- entropy_all + entropy_type 
          }
        }
        return(entropy_all)
      }
    }	
  }else{
    print("Cell type not found!")
    return(NA)
  }
}
