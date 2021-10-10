#' entropy
#'
#' @description Produces a dataframe of entropies for all reference cells or an entropy number for the whole image if a radius is not supplied. 
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param types_of_interest Cell.Types of interest, the first celltype is considered as referece celltype
#' @param radius (OPTIONAL) The maximum radius around a reference cell for another cell to be considered an interaction.
#' @param column String specifying the column the markers are from
#' @importFrom SummarizedExperiment colData 
#' @return A dataframe or a number depending on the argument radius
#' @export
entropy <- function(sce_object, types_of_interest, column, radius = NULL){
  if((types_of_interest[1] %in% colData(sce_object)[,column]) && 
     any(types_of_interest[-1] %in% colData(sce_object)[,column])){
    if (!is.null(radius)){
      reference_marker <- types_of_interest[1]
      target_marker<-types_of_interest
      n_cells <- number_of_cells_within_radius(sce_object, reference_marker = reference_marker,
                                               target_marker = target_marker,radius = radius, column = column)
      n_cells.df <- n_cells[[reference_marker]]
      n_cells.df[,"total"] <- 0
      
      for (i in target_marker){
        n_cells.df[,paste(i,"_log10",sep = "")] <- log10(n_cells.df[,i])
        n_cells.df[,"total"] <- n_cells.df[,"total"]+n_cells.df[,i]
      }
      n_cells.df[,"total_log"] <- log10(n_cells.df[,"total"])
      n_cells.df[n_cells.df == -Inf] <- 0
      
      for (i in target_marker){
        n_cells.df[which(n_cells.df[,"total"] != 0), paste(i,"ratio",sep = "")] <- n_cells.df[,i]/n_cells.df[,"total"]
        n_cells.df[which(n_cells.df[,"total"] == 0), paste(i,"ratio",sep = "")] <- 0
        n_cells.df[,paste(i,"_entropy",sep = "")] <- n_cells.df[, paste(i,"ratio",sep = "")] *(n_cells.df[,paste(i,"_log10",sep = "")] - n_cells.df[,"total_log"])
      }
      n_cells.df[,"entropy"] <- -(rowSums(n_cells.df[,grepl("_entropy",colnames(n_cells.df))]))  

      return(n_cells.df)
    } else{
      entropy_all <- 0
      data <- data.frame(colData(sce_object))
      data <- data[which(data[,column] %in% types_of_interest),]
      n_all <- dim(data)[1]
      
      if (n_all == 0){
        return(0.0)
      } else{
        for (type in types_of_interest){
          type_data <- data[which(data[,column]==type),]
          n_type <- dim(type_data)[1]
          if (n_type != 0){
            p_type <- n_type/n_all
            entropy_type <- -(p_type)*log10(p_type)
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
