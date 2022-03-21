#' The ratio of tumour border cell count and tumour cell count
#'
#' @description Calculates the ratio of the tumour border cell count and the
#'   total tumour cell count in an image. The tumour border cells are detected
#'   by the default `identify_bordering_cells` function. If the ratio is high,
#'   it means that most tumour cells are identified as bordering cells. This
#'   means there is no clear tumour clusters.
#' @param sce_object SingleCellExperiment object in the form of the output of
#'   format_image_to_sce.
#' @param cell_type_of_interest String. The cell type that the user wants to
#'   determine a cluster of.
#' @param feature_colname String. The column that contains the
#'   `cell_type_on_interest`.
#'
#' @return A number is returned.
#' @export
#' @examples 
#' R_BT(SPIAT::defined_image, cell_type_of_interest = "Tumour", "Cell.Type")
#' 
R_BT <- function(sce_object, cell_type_of_interest, feature_colname){
  
  # identify the bordering cells
  sce_border <- identify_bordering_cells(sce_object, reference_cell = cell_type_of_interest, 
                                         feature_colname = feature_colname, 
                                         ahull_alpha = 40,
                                         n_to_exclude = 0)
  
  # count the number of bordering cells and tumour cells
  n_tumour <- count_category(sce_border, cat = cell_type_of_interest,
                             feature_colname = feature_colname)
  n_border <- count_category(sce_border, "Border","Region")
  # calculate the ratio
  r <- n_border/n_tumour
  
  return(r)
}
