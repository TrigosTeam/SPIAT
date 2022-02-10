#' R_BT
#'
#' @param sce_object 
#' @param cell_types_of_interest 
#' @param feature_colname 
#'
#' @return
#' @export
R_BT <- function(sce_object, cell_types_of_interest, feature_colname){
  
  # define a function to get the number of certain name under a certain column
  count_category <- function(sce_object, 
                             cat = cell_types_of_interest, 
                             feature_colname = feature_colname){
    data <- data.frame(colData(sce_object))
    count_table <- table(data[feature_colname])
    count <- unname(count_table[match(cat, names(count_table))])
    return(count)
  }
  
  # identify the bordering cells
  sce_border <- identify_bordering_cells(sce_object, reference_cell = cell_types_of_interest, 
                                         feature_colname = feature_colname, 
                                         ahull_alpha = 40,
                                         n_to_exclude = 0)
  
  # count the number of bordering cells and tumour cells
  n_tumour <- count_category(sce_border, cat = cell_types_of_interest,
                             feature_colname = feature_colname)
  n_border <- count_category(sce_border, "Border","Region")
  # calculate the ratio
  r <- n_border/n_tumour
  
  return(r)
}


