#' AUC_of_cross_function
#'
#' @description Calculate the difference of area under curve between two curves 

#' @param df.cross Dataframe containing the positions of the two curves 
#' @importFrom pracma trapz
#' @export


AUC_of_cross_function <- function(df.cross){
  AUC <- pracma::trapz(df.cross$r ,df.cross$border) - pracma::trapz(df.cross$r ,df.cross$theo)
  
  # get the cross k result image size
  X <- max(df.cross$theo)
  Y <- max(df.cross$r)
  
  # calculate normalised AUC
  n_AUC <- AUC/(X*Y)
  return(n_AUC)
}
