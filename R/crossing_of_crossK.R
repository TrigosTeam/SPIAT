#' crossing_of_crossK
#'
#' @description Determine if there is a crossing in the cross K curves, to
#'   further detect the existence of potential immune rings. --requires
#'   additional explanation/reference. Cite own paper?
#'
#' @param df.cross Data.frame. The output of
#'   \code{\link{calculate_cross_functions}}. Containing the positions of the
#'   two curves. Columns contain "r", "border" and "theo".
#'
#' @return A number. The percentage of the crossing position of the specified
#'   distance.
#' @export
#'
#' @examples
#' df_cross <- calculate_cross_functions(SPIAT::defined_image, method="Kcross",
#'               cell_types_of_interest = c("Tumour","Immune3"),
#'               feature_colname ="Cell.Type", dist = 100)
#' crossing_of_crossK(df_cross)

crossing_of_crossK <- function(df.cross){
    df.cross$sign <- df.cross$theo - df.cross$border
    change_of_sign <- diff(sign(df.cross$sign[-1]))
    ix <- which(change_of_sign != 0)
    n <- dim(df.cross)[1]
    if (length(ix) == 1 && ix/n > 0.04){
        methods::show("Crossing of cross K function is detected for this image, indicating a potential immune ring.")
        perc <- round(ix/n *100,2)
        methods::show(paste("The crossing happens at the", 
                   perc, "% of the specified distance.", sep = ""))
    }
    else ix <- NA
    return(round(1-ix/n,2))
}