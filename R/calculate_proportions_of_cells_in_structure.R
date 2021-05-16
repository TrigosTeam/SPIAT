#' calculate_proportions_of_cells_in_structure
#'
#' @description Calculate the proprtion of immune cells in each defined tumour structure, return a dataframe 
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @export

calculate_proportions_of_cells_in_structure <- function(sce_object){
  
  data <- data.frame(colData(sce_object))
  
  # proportion of infiltrated immune cells that are in the tumour region
  p_i_in <- length(which(data$Structure == "Infiltration"))/length(which(data$Region == "In"))
  # proportion of immune cells that are in the tumour side of the invasive front
  p_i_i.f.in <- length(which(data$Structure == "I.f.immune.in"))/
    length(which(data$Structure == "Invasive.front.in"))
  # proportion of immune cells of the invasive front that are outside of the tumour 
  p_i_i.f.out <- length(which(data$Structure == "I.f.immune.out"))/
    length(which(data$Structure == "Invasive.front.out"))
  # proportion of exclusice immune cells that are out of the tumour region
  p_i_out <- length(which(data$Structure == "Exclusion"))/length(which(data$Region == "Out"))
  
  #save in dataframe
  Infiltration <- p_i_in
  Exclusion <- p_i_out
  Invasive.front.in <- p_i_i.f.in
  Invasive.front.out <- p_i_i.f.out
  proportion.data <- data.frame(Infiltration, Exclusion, Invasive.front.in, Invasive.front.out)
  
  return(proportion.data)
}

