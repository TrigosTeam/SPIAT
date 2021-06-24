#' calculate_proportions_of_cells_in_structure
#'
#' @description Calculate the proprtion of immune cells in each defined tumour structure, return a dataframe 
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @export

calculate_proportions_of_cells_in_structure <- function(sce_object){
  
  data <- data.frame(colData(sce_object))
  
  # proportion of infiltrated immune cells that are in the tumour region
  p_i_in <- length(which(data$Structure == "Infiltrated.immune"))/length(which(data$Region == "Inside"))
  # proportion of immune cells that are in the tumour side of the invasive front
  p_i_i.f.in <- length(which(data$Structure == "Internal.margin.immune"))/
    length(which(data$Structure == "Internal.margin"))
  # proportion of immune cells of the invasive front that are outside of the tumour 
  p_i_i.f.out <- length(which(data$Structure == "External.margin.immune"))/
    length(which(data$Structure == "External.margin"))
  # proportion of exclusice immune cells that are out of the tumour region
  p_i_out <- length(which(data$Structure == "Stromal.immune"))/length(which(data$Region == "Outside"))
  
  #save in dataframe
  Infiltrated.immune <- p_i_in
  Stromal.immune <- p_i_out
  Internal.margin.immune <- p_i_i.f.in
  External.margin.immune <- p_i_i.f.out
  proportion.data <- data.frame(Infiltrated.immune, Stromal.immune, 
                                Internal.margin.immune, External.margin.immune)
  
  return(proportion.data)
}

