#' Find and plot cluster model parameters of multiple images in a boxplot
#' 
#' fit_cluster_model_multiple_images takes in a list containing image point patterns, specific phenotypes
#' may be specified. The function outputs a list containing vectors containing kappa, estimated point intensity and mean
#' cluster size. The function by default displays the data as boxplots
#' 
#' Inputs - 
#' image_point_patterns: list of ppp objects to be analysed
#' plot: Parameter to indicate whether plotting of data is necessary 
#' phenotypes_of_interest: Cell phenotypes to be considered when fitting the cluster model
#' 
#' @param image_point_patterns list containing point pattern objects to be compared
#' @param plot variable to indicate whether to plot the cluster model data, default is to plot, set to 
#' NULL if plot is not required
#' @param phenotypes_of_interest string vector containing name of cell phenotype of interest, can be a single or multiple phenotypes. 
#' By default it is not specified 
#' 

fit_cluster_model_multiple_images <- function(image_point_patterns = NULL, phenotypes_of_interest = NULL, plot = 1){
  # Determine number of images to iterate through 
  i <- 1: length(image_point_patterns)
  
  # Declare empty variables to be filled with data
  empty_vector <- vector(mode = "numeric", length = length(image_point_patterns))
  summary_values <- list(kappa = empty_vector, mean_cluster_size = empty_vector, estimated_intensity = empty_vector)
  
  # Subset point patterns if necessary 
  if(!is.null(phenotypes_of_interest))
  {
    temp_point_patterns <- vector(mode = "list", length = length(image_point_patterns))
    for(val in seq_along(i))
    {
      temp_point_patterns[[val]] <- find_subset(point_pattern = image_point_patterns[[val]], cell_phenotypes_of_interest = phenotypes_of_interest)
    }
    image_point_patterns <- temp_point_patterns
  }
  
  # Fit cluster model and extract parameters for each image into vectors
  for (val in seq_along(i))
  {
    fitted_model <- fit_cluster_model(image_point_patterns[[val]])
    summary_values$kappa[val] <- fitted_model$kappa
    summary_values$mean_cluster_size[val] <- fitted_model$mean_cluster_size
    summary_values$estimated_intensity[val] <- fitted_model$estimated_intensity
  }
  
  # Plot boxplots of all data if required
  if(!is.null(plot))
  {
    boxplot(summary_values$kappa, main = "Summary of kappa values", xlab = "kappa", horizontal = TRUE)
    boxplot(summary_values$mean_cluster_size, main = "Summary of cluster size", xlab = "Mean Cluster Size", horizontal = TRUE)
    boxplot(summary_values$estimated_intensity, main = "Summary of point intensity", xlab = "Estimated intensity", horizontal = TRUE)
  }
  
  return(summary_values)
}
