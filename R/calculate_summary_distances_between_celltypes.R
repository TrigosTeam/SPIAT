#'calculate_summary_distances_between_celltypes
#'
#'@description Returns the mean, median and standard deviation of the
#'  minimum/pairwise distances between phenotypes.
#'
#'@param df Data.frame containing the distance output between cell types. The
#'  functions that generate the distances can be
#'  \code{\link{calculate_minimum_distances_between_celltypes}} and
#'  \code{\link{calculate_pairwise_distances_between_celltypes}}.
#' @import dplyr
#' @return A data frame is returned
#' @examples
#' # for pairwise dist
#' pairwise_dist <- calculate_pairwise_distances_between_celltypes(
#' SPIAT::defined_image, cell_types_of_interest = c("Tumour","Immune1"), 
#' feature_colname = "Cell.Type")
#' summary_distances <- calculate_summary_distances_between_celltypes(pairwise_dist)
#' 
#' # for minimum dist
#' min_dists <- calculate_minimum_distances_between_celltypes(
#' SPIAT::defined_image, cell_types_of_interest = c("Tumour","Immune1"), 
#' feature_colname = "Cell.Type")
#' summary_distances <- calculate_summary_distances_between_celltypes(min_dists)
#'@export

calculate_summary_distances_between_celltypes <- function(df) {
    Pair <- Distance <- NULL
    # summarise the results
    summarised_dists <- df %>% 
        group_by(Pair) %>%
        summarise(mean(Distance, na.rm = TRUE), min(Distance, na.rm = TRUE), 
                  max(Distance, na.rm = TRUE),
                  stats::median(Distance, na.rm = TRUE), 
                  stats::sd(Distance, na.rm = TRUE))
    
    summarised_dists <- data.frame(summarised_dists)
    colnames(summarised_dists) <- c("Pair" , "Mean", "Min", "Max", 
                                    "Median", "Std.Dev")
    
    for (i in seq(1,dim(summarised_dists)[1])){
        cellnames <- strsplit(summarised_dists[i,"Pair"],"/")[[1]]
        summarised_dists[i, "Reference"] <- cellnames[1]
        summarised_dists[i, "Target"] <- cellnames[2]
    }
    return(summarised_dists)
}
