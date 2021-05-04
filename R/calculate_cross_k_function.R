#' Perform cross-K function
#' 
#' calculate_cross_k_function, Calculates Cross-K function between 2 specified cell phenotype groups.
#' The phenotype groups can contain multiple or single phenotypes. Options are given to name the two groups and to plot the function.
#' 
#' Inputs - 
#' point_pattern: marked ppp object representing intial image
#' phenotype_group1: String containing the names of the phenotypes to be considered, can be singular or multiple.
#' phenotype_group2: See phenotype_group_1
#' plot: Indicates whether to plot Cross-K function by checking for numeric value. By default set to plot graph (plot = 1), if no plot required set to NULL
#' names: String vector containing desired names for the two groups of phenotypes, largely for plotting purposes. If not specified groups are named Group1 and Group2.
#' 
#' @description Calculate Cross-K function for 2 groups of cell phenotypes of interest
#' 
#' @param point_pattern ppp object representing image
#' @param phenotype_group1 Cell phenotypes to be compared from, can be a single phenotype or multiple
#' @param phenotype_group2 Cell phenotypes to be compared to, can be a single phenotype or multiple
#' @param plot Dictate whether to plot Cross-K function 
#' @param names String vector containing the names of the phenotype groups
#' @importFrom spatstat.core ppm Kcross.inhom 
#' @importFrom spatstat.geom superimpose

calculate_cross_k_function <- function(point_pattern, phenotype_group1, phenotype_group2, names = NULL, plot = 1) {
  # Create subsets for Cross-K function comparison
  subset1 <- find_subset(point_pattern, cell_phenotypes_of_interest = phenotype_group1)
  subset2 <- find_subset(point_pattern, cell_phenotypes_of_interest = phenotype_group2)
  
  # Assign relevant marker names
  if(!is.null(names))  {
    marker1 <- names[1]
    marker2 <- names[2]
  }
  # Use supplied marker names
  else  {
    marker1 <- "Group 1"
    marker2 <- "Group 2"
  }
  
  # Assign markers to subsets
  marks(subset1) <- as.factor(marker1)
  marks(subset2) <- as.factor(marker2)
  
  # Create model based on point subsets to be used when calculating Cross-K
  group1 <- ppm(subset1)
  group2 <- ppm(subset2)
  
  # Create new point pattern with chosen points
  new_point_pattern <- superimpose(subset1, subset2)
  
  # Calculate Cross-K
  cross_k_function <- Kcross.inhom(new_point_pattern, marker1, marker2, group1, group2, correction = c("border"))
  
  # Plot if required
  if(!is.null(plot)){
    if(plot == 2){
      plot(cross_k_function, main = NULL, xlab = "Radius")
    }
    else {
      plot(cross_k_function, main = "Cross-K function", xlab = "Radius")
    }
  }
  
  # Return Cross-K function data
  return(cross_k_function)
}
  
