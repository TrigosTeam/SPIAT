#' crossing_of_crossK
#'
#' @description Determine if there is a crossing in the cross K curves, to further 
#' detect the existance of potential immune ring.
#' 
#' @param df.cross Dataframe containing the positions of the two curves 
#'
#' @return A numeric The percentage of the crossing position of the specified distance.
#' @export
#'
#' @examples
crossing_of_crossK <- function(df.cross){
  df.cross$sign <- df.cross$theo - df.cross$border
  change_of_sign <- diff(sign(df.cross$sign[-1]))
  ix <- which(change_of_sign != 0)
  n <- dim(df.cross)[1]
  if (length(ix) == 1 && ix/n > 0.04){
    print("Crossing of cross K function is detected for this image, indicating a potential immune ring.")
    perc <- round(ix/n *100,2)
    print(paste("The crossing happens at the", perc, "% of the specified distance."), sep = "")
  }
  else ix <- NA
  return(round(1-ix/n,2))
}
