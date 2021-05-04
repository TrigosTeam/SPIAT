#' Determine point subset data 
#' 
#' find_subset returns unmarked point pattern representing the phenotype/phenotypes 
#' of interest from the initial image
#' 
#' Inputs - 
#' point_pattern: Marked point pattern representing data
#' cell_phenotypes_of_interest: String representation of the names of the phenotypes of interest
#' 
#' Outputs - 
#' subset_of_interest: unmarked ppp object which is subset of points representing the phenotype specified, 
#' 
#' @param point_pattern ppp object representing initial image
#' @param cell_phenotypes_of_interest Phenotype/Phenotypes of interest
#' @importFrom stringr str_detect 
#' @importFrom spatstat.geom superimpose
find_subset <- function(point_pattern, cell_phenotypes_of_interest){
  # Split initial ppp object into subsets based on marks
  subset <- split(point_pattern)
  
  # Check if there are multiple phenotypes of interest
  if(length(cell_phenotypes_of_interest) == 1)
  {
    # Check whether the phenotype has a comma to indicate multiple markers and act accordingly
    if(stringr::str_detect(cell_phenotypes_of_interest, ",")){
      subset_name <- paste("subset$`", cell_phenotypes_of_interest, "`", sep = "")
    }
    else{
      subset_name <- paste("subset$", cell_phenotypes_of_interest, sep = "")
    }
    subset_of_interest <- eval(parse(text = subset_name))
  }
  # Check for multiple phenotypes of interest and act accordingly to superimpose the subsets required
  else if(length(cell_phenotypes_of_interest) > 1)
  {
    for (j in 1:length(cell_phenotypes_of_interest))
    {
      if(stringr::str_detect(cell_phenotypes_of_interest[j], ",")){
        subset_name <- paste("subset$`", cell_phenotypes_of_interest[j], "`", sep = "")
      }
      else{
        subset_name <- paste("subset$", cell_phenotypes_of_interest[j], sep = "")
      }
      subset_of_interest <- superimpose(subset_of_interest, eval(parse(text = subset_name)))
    }
  }
  return(subset_of_interest)
}

