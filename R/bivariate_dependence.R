#' bivariate_dependence
#'
#' @description find the dependence between two types of cells
#' 
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param method String that is the method for dependence calculation
#' @param phenotypes Vector of phenotypes of interest
#' @param merge Vector of phenotypes to be merged
#' @param merge_name String that is the name of the merged phenotype
#' @import SingleCellExperiment
#' @export 

# colData() is in package 'SummarizedExperiment' but imported by SingleCellExperiment

bivariate_dependence <- function(sce_object, method = "Kcross", 
                                 phenotypes, merge = NULL, 
                                 merge_name = NULL) {
  
  formatted_data <- colData(sce_object)

  #CHECK
  if (!all(phenotypes %in% formatted_data$Phenotype)) {
    stop("phenotype not found")
  }
  
  # CHECK if there are any columns to be merged
  if (is.null(merge) == FALSE){
    formatted_data[which(formatted_data$Phenotype %in% merge),]$Phenotype <- merge_name
    # update phenotypes of interest
    phenotypes <- c(setdiff(phenotypes, merge), merge_name)
  }
  
  # get x, y coordinates and phenotypes from formatted_data
  x <- formatted_data$Cell.X.Position
  y <- formatted_data$Cell.Y.Position
  marks <- formatted_data$Phenotype
  
  # get windows
  x_summary <- summary(formatted_data$Cell.X.Position)
  x_min <- as.numeric(x_summary[1])
  x_max <- as.numeric(x_summary[6])
  y_summary <- summary(formatted_data$Cell.Y.Position)
  y_min <- as.numeric(y_summary[1])
  y_max <- as.numeric(y_summary[6])
  
  
  # format sce to ppp
  ppp_object <- ppp(x, y, c(x_min,x_max),c(y_min,y_max), marks = marks)
  ppp_object$marks <- as.factor(ppp_object$marks)
  
  if (method == "Gcross"){
    plot(Gcross(ppp_object, phenotypes[1],phenotypes[2]), main = "cross G function")
  }
  if (method == "Kcross"){
    plot(Kcross(ppp_object, phenotypes[1],phenotypes[2]), main = "cross K function")
  }
  if (method == "Lcross"){
    plot(Lcross(ppp_object, phenotypes[1],phenotypes[2]), main = "cross L function")
  }
  if (method == "Jcross"){
    plot(Jcross(ppp_object, phenotypes[1],phenotypes[2]), main = "cross J function")
  }
}

