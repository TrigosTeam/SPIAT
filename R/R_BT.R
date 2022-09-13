#' The ratio of border cell count to cluster cell count
#'
#' @description Calculates the ratio of the bordering cell count and the
#'   total to-be-clustered cell count in an image. The bordering cells are detected
#'   by the default \code{\link{identify_bordering_cells}} function. If the
#'   ratio is high, it means that most cells to be clustered are identified as bordering
#'   cells. This means there is no clear clusters.
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param cell_type_of_interest String. The cell type that the user wants to
#'   determine a cluster of.
#' @param feature_colname String. The column that contains the cell type to be clustered.
#'
#' @return A number is returned.
#' @export
#' @examples
#' R_BC(SPIAT::defined_image, cell_type_of_interest = "Tumour", "Cell.Type")
#' 
R_BC <- function(spe_object, cell_type_of_interest, feature_colname){
    # identify the bordering cells
    spe_border <- identify_bordering_cells(spe_object, 
                                           reference_cell=cell_type_of_interest,
                                           feature_colname = feature_colname, 
                                           ahull_alpha = 40,
                                           n_to_exclude = 0)
    
    # count the number of bordering cells and tumour cells
    n_cluster_cell <- count_category(spe_border, cat = cell_type_of_interest,
                               feature_colname = feature_colname)
    n_border_cell <- count_category(spe_border, "Border","Region")
    # calculate the ratio
    r <- n_border_cell/n_cluster_cell
    
    return(r)
}
