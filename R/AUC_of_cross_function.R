#' AUC_of_cross_function
#'
#' @description Calculate the difference of area under curve between two curves 
#'
#' @param df.cross Dataframe containing the positions of the two curves 
#' @import pracma
#' @export

#sce_object = formatted_defined
#reference_cell = "MEL"
#n_of_polygons = 2
#buffer_width = 30
#ahull_alpha = 60

# colData() is in package 'SummarizedExperiment' but imported by SingleCellExperiment

AUC_of_cross_function <- function(df.cross){
  AUC <- pracma::trapz(df.cross$r ,df.cross$border) - pracma::trapz(df.cross$r ,df.cross$theo)
  return(AUC)
}