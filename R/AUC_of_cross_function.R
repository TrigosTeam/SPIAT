#' AUC_of_cross_function
#'
#' @description Calculate the difference of area under curve between two curves 
#'
#' @param sce_object SingleCellExperiment object that was used for cross K calculation
#' @param df.cross Dataframe containing the positions of the two curves 
#' @importFrom pracma trapz
#' @export


AUC_of_cross_function <- function(sce_object, df.cross){
  AUC <- pracma::trapz(df.cross$r ,df.cross$border) - pracma::trapz(df.cross$r ,df.cross$theo)
  
  # get the image size
  data <- data.frame(colData(sce_object))
  X <- max(data$Cell.X.Position)
  Y <- max(data$Cell.Y.Position)
  
  # calculate normalised AUC
  n_AUC <- AUC/(X*Y)
  return(n_AUC)
}
