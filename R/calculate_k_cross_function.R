#' Perform cross-K function
#' 
#' @description calculate cross-K function for 2 chosen phenotypes of interest
#' 
#' @param point_pattern ppp object representing image
#' @param cell_phenotypes_of_interest cell phenotypes to be compared
#' 
#' @importFrom spatstat.core Kcross.inhom
#' @importFrom spatstat.geom superimpose
#' @importFrom stringr str_detect

calculate_k_cross_function <- function(point_pattern, cell_phenotypes_of_interest = NULL)
{
  
  subset <- split(point_pattern)
  marker1 <- cell_phenotypes_of_interest[1]
  marker2 <- cell_phenotypes_of_interest[2]
  
  if(stringr::str_detect(marker1, ",")){
    subset_name1 <- paste("subset$`", marker1, "`", sep = "")
  }
  else{
    subset_name1 <- paste("subset$", marker1, sep = "")
  }
  
  if(stringr::str_detect(marker2, ",")){
    subset_name2 <- paste("subset$`", marker2, "`", sep = "")
  }
  else{
    subset_name2 <- paste("subset$", marker2, sep = "")
  }
  
  subset1 <- eval(parse(text = subset_name1))
  subset2 <- eval(parse(text = subset_name2))
  
  marks(subset1) <- as.factor(marker1)
  marks(subset2) <- as.factor(marker2)
  
  group1 <- ppm(subset1)
  group2 <- ppm(subset2)
  
  new_point_pattern <- superimpose(subset1, subset2)
  
  k_cross_function <- Kcross.inhom(new_point_pattern, marker1, marker2, group1, group2)
  plot(k_cross_function)
  
  return(k_cross_function)
}
  
