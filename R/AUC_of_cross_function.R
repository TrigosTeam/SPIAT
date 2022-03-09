#' The difference in AUC of the cross function curves
#'
#' @description Calculate the difference of area under the curve (AUC) between
#'   two curves, normalised by the total area of the graph.

#' @param df.cross Data.Frame. The output the the cross functions. Containing
#'   the positions of the two curves.
#' @export


AUC_of_cross_function <- function(df.cross){
  AUC <- pracma::trapz(df.cross$r ,df.cross$border) - 
    pracma::trapz(df.cross$r,df.cross$theo)
  
  # get the cross k result image size
  X <- max(df.cross$theo)
  Y <- max(df.cross$r)
  
  # calculate normalised AUC
  n_AUC <- AUC/(X*Y)
  return(n_AUC)
}
