#' calculate_proportions_of_immune_cells_in_structure
#'
#' @description Calculate the proprtion of immune cells in each defined tumour structure, return a data_localframe 
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param immune_cells Vector of immune cells to consider
#' @param column Column to extract cell types from
#' @export


calculate_proportions_of_immune_cells_in_structure <- function(sce_object, immune_cells, column){
  
  data_local <- data.frame(colData(sce_object))
  
  ##Relative to all cells
  proportions <- vector()
  for(immune_cell_type in immune_cells){
    temp <- data_local[data_local[,column] == immune_cell_type,]
    p_i_in <- length(which(temp$Structure == "Infiltrated.immune"))/length(which(data_local$Region == "Inside"))
    p_i_i.f.in <- length(which(temp$Structure == "Internal.margin.immune"))/
      length(which(data_local$Structure == "Internal.margin"))
    p_i_i.f.out <- length(which(temp$Structure == "External.margin.immune"))/
      length(which(data_local$Structure == "External.margin"))
    p_i_out <- length(which(temp$Structure == "Stromal.immune"))/length(which(data_local$Region == "Outside"))
    
    proportions <- rbind(proportions,
                         c(immune_cell_type, "To_all_cells", p_i_in, p_i_i.f.in, p_i_i.f.out, p_i_out))
  }

  proportions <- as.data.frame(proportions)

  ##Relative to all immune cells in each area
  data_local_immune <- data_local[data_local[,column] %in% immune_cells,]
  for(immune_cell_type in immune_cells){
    temp <- data_local[data_local[,column] == immune_cell_type,]
    p_i_in <- length(which(temp$Structure == "Infiltrated.immune"))/length(which(data_local_immune$Structure == "Infiltrated.immune"))
    p_i_i.f.in <- length(which(temp$Structure == "Internal.margin.immune"))/
      length(which(data_local_immune$Structure == "Internal.margin.immune"))
    p_i_i.f.out <- length(which(temp$Structure == "External.margin.immune"))/
      length(which(data_local_immune$Structure == "External.margin.immune"))
    p_i_out <- length(which(temp$Structure == "Stromal.immune"))/length(which(data_local_immune$Structure == "Stromal.immune"))
    
    proportions <- rbind(proportions,
                         c(immune_cell_type, "To_immune_cells", p_i_in, p_i_i.f.in, p_i_i.f.out, p_i_out))
  }
  colnames(proportions) <- c("Cell.Type2", "Relative_to", "P.Infiltrated.Immune",
                             "P.Internal.Margin.Immune",
                             "P.External.Margin.Immune", "P.Stromal.Immune")
  proportions$P.Infiltrated.Immune <- as.numeric(proportions$P.Infiltrated.Immune)
  proportions$P.Stromal.Immune <- as.numeric(proportions$P.Stromal.Immune)
  proportions$P.Internal.Margin.Immune <- as.numeric(proportions$P.Internal.Margin.Immune)
  proportions$P.External.Margin.Immune <- as.numeric(proportions$P.External.Margin.Immune)
  
  return(proportions)
}
