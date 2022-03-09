#' check for initial underlying distribution (cluster or Homogenous Poisson)
#' 
#' check_distribution_types Takes in point pattern object and performs a quadrat test and clark evans test
#' to determine the underlying pattern. The results and p value will be printed to console
#' 
#' Input - 
#' point_pattern: ppp object representing image to be analysed
#' quadrat_dim: Numeric vector containing the dimensions for image segmentation when performing the quadrat test. 
#' default is 3x3 segementation in x,y direction
#' phenotypes_of_interest: String vector containing the names of the cell phenotypes 
#' to be examined in the case the entire point pattern is not relevant
#' 
#' Output - 
#' test_results: List containing 2 htest objects containing test data 
#' 
#' @param point_pattern sample point pattern to be analysed
#' @param quadrat_dim vector containing number of quadrants in x,y direction respectively
#' @param phenotypes_of_interest String describing the phenotypes to be examined 
#' @importFrom spatstat.core quadrat.test clarkevans.test

check_distribution_type <- function(point_pattern, quadrat_dim = c(3,3), phenotypes_of_interest = NULL){
  # Find point pattern subset if required
  if(is.null(phenotypes_of_interest)){
    point_pattern <- find_subset(point_pattern, phenotypes_of_interest)
  }
  
  # Perform Quadrat Test 
  quadrat_test_results <- quadrat.test(point_pattern, quadrat_dim[1],
                                       quadrat_dim[2], alternative = "two.sided")
  # Print quadrat test results 
  if(quadrat_test_results$p.value < 0.05){
    print(paste("Reject null distribution (Complete Spatial Randomness), p-value:", signif(quadrat_test_results$p.value, 4)))
    print("Suggests Inhomogenous Poisson Distribution")
  }
  else if(quadrat_test_results$p.value >= 0.05){
    print(paste("Accept null distribution (Complete Spatial Randomness), p-value:", signif(quadrat_test_results$p.value, 4)))
    print("Suggests Homogenous Poisson Distribution")
  }
  
  # Perform Clark Evans test 
  clark_evans_results <- clarkevans.test(point_pattern, alternative="clustered", correction = "none")
  # Print Clark Evans Test results
   if (clark_evans_results$statistic < 1){
    print(paste("Reject null distribution (Complete Spatial Randomness), p-value:", signif(clark_evans_results$p.value, 4)))
    print("Suggests clustered point pattern")
  }
  else if (clark_evans_results$statistic >= 1){
    print(paste("Accept null distribution (Complete Spatial Randomness), p-value:", signif(clark_evans_results$p.value, 4)))
    print("Suggests ordered point pattern")
  }
  
  # Create and return list with test results
  test_results <- list("quadrat" = quadrat_test_results, "clarkevans" = clark_evans_results)
  return(test_results)
}
