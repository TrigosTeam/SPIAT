#' format_colData_to_sce
#'
#' @description Formats a dataframe of colData into a singlecellexperiment class 
#' where the count assay is empty for
#' every cell (columns), and cell phenotype, x and y coordinates are stored under colData
#' for the purpose of passing dataframe into a function requiring sce_object
#'
#' @export
#' @param colData Dataframe that will be the colData of the sce object
#' @import SingleCellExperiment


format_colData_to_sce <- function(colData) {
  colData[,"pseudo"] <- 0
  assay_data <- colData[,"pseudo"]
  assay_rownames <- "pseudo"
  assay_colnames <- rownames(colData)
  
  #transpose the matrix so every column is a cell and every row is a marker
  assay_data_matrix <- as.matrix(assay_data)
  colnames(assay_data_matrix) <- NULL
  rownames(assay_data_matrix) <- NULL
  assay_data_matrix_t <- t(assay_data_matrix)
  
  sce <- SingleCellExperiment(assays = list(counts = assay_data_matrix_t))
  
  rownames(sce) <- assay_rownames
  colnames(sce) <- assay_colnames
  
  #Assign the phenotype, X and Y positions as the colData
  coldata_phenotype <- colData[,"Phenotype"]
  coldata_Xpos <- colData[,"Cell.X.Position"]
  coldata_Ypos <- colData[,"Cell.Y.Position"]
  colData(sce)$Phenotype <- coldata_phenotype
  colData(sce)$Cell.X.Position <- coldata_Xpos
  colData(sce)$Cell.Y.Position <- coldata_Ypos
  return(sce)
}