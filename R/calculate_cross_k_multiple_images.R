#' Calculate K-cross function for multiple images and plot
#' 
#' Inputs -
#' image_point_patterns: List containing point pattern objects of the images of interest
#' phenotype_group1: Names of phenotypes of interest for group 1
#' phenotype_group2: Names of phenotypes of interest for group 2
#' names: String vector containing desired names for the two groups of phenotypes, largely for plotting purposes. If not specified groups are named Group1 and Group2.
#' plot: Dictates whether to plot cross K-functions calculated, set to NULL if plotting is not required
#' 
#' Outputs - 
#' k_cross_functions: list containing the data required to plot multiple cross-K functions
#' 
#' @param image_point_patterns Point patterns to be investigated
#' @param phenotype_group1 Names of phenotypes of interest (Group1)
#' @param phenotype_group2 Names of phenotypes of interest (Group2)
#' @param plot Dictate whether plotting of function is required 
#' @param names String vector containing the names of the phenotype groups
#' 
#' @importFrom graphics title 
calculate_cross_k_multiple_images <- function(image_point_patterns, phenotype_group1, phenotype_group2, names = NULL, plot = 2){
  i <- 1: length(image_point_patterns)
  cross_k_functions <- vector(mode = "list", length(image_point_patterns))
  # Iterate through images and calculate cross k function to be stored in k_cross_function
  for(val in i){
    cross_k_functions[[val]] <- calculate_cross_k_function(image_point_patterns[[val]], phenotype_group1, phenotype_group2, plot = plot)
    title(main = paste("Cross K-Function - Point Pattern", val))
  }
  return(cross_k_functions)
}
