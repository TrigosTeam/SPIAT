#' calculate_k_function 
#' 
#' @description Returns K function of cell distribution (inhomogenous distributions)
#' 
#' @param point_pattern ppp object used as input data for K function 
#' @param determine_confidence indicate the confidence interval of the K function calculated
#' @param cell_phenotype_of_interest string supplying the phenotype of interest

calculate_k_function <- function(point_pattern, cell_phenotype_of_interest = NULL, determine_confidence = NULL, plot_envelope = NULL){
 library(stringr)
  
  if(!is.null(cell_phenotype_of_interest)) {
    subset <- split(point_pattern)
    if(str_detect(cell_phenotype_of_interest, ",")){
      subset_name <- paste("subset$`", cell_phenotype_of_interest, "`", sep = "")
    }
    else{
      subset_name <- paste("subset$", cell_phenotype_of_interest, sep = "")
    }
    subset_of_interest <- eval(parse(text = subset_name))
    point_pattern <- subset_of_interest
  } 
  D <- density(point_pattern)
  k_function <- Kinhom(point_pattern, correction = "border", D)
  
  
  if (!is.null(determine_confidence)) {
    confidence <- lohboot(point_pattern, "Kinhom")
    plot(confidence)
  }
  if (!is.null(plot_envelope)) {
    envelope_k <- envelope(point_pattern, Kinhom, funargs = c(D, correction = "border"),  simulate = expression(rpoispp(density(subset2))), nsim= 49)
   
    plot(envelope_k, main = "K-function")
  } 
  else {
    plot(k_function, main = "K-function")
  }
  return(k_function)
}
