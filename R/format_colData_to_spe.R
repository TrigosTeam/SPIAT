#' format_colData_to_spe
#'
#' @description Format a data frame into a SpatialExperiment class where the
#'   count assay is empty every cell (columns), cell phenotypes are stored under
#'   colData() and cell coordinates are stored under spatialCoords().
#'
#' @param df Data frame that contains cell coordinates, phenotypes (if
#'   available) and other cell properties. The rownames should be cell ID
#' @return An SpatialExperiment object
#' @examples
#' df <- data.frame(row.names = c("Cell_1", "Cell_2"), Cell.X.Position = c(2,5),
#' Cell.Y.Position = c(3.3, 8), Phenotypes = c("CD3", "CD3,CD8"))
#' spe <- format_colData_to_spe(df)
#' @export

format_colData_to_spe <- function(df) {
    
    #CHECK
    if (dim(df)[1]==0){
        stop("No data in the dataframe")
    } 
    
    assay_data <- rep(0, dim(df)[1])
    assay_rownames <- "pseudo"
    assay_colnames <- rownames(df)
    
    #transpose the matrix so every column is a cell and every row is a marker
    assay_data_matrix <- as.matrix(assay_data)
    colnames(assay_data_matrix) <- NULL
    rownames(assay_data_matrix) <- NULL
    assay_data_matrix_t <- t(assay_data_matrix)
    
    spe <- SpatialExperiment::SpatialExperiment(
        assays = assay_data_matrix_t, colData = df, 
        spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position"))
    
    rownames(spe) <- assay_rownames
    colnames(spe) <- assay_colnames
    
    return(spe)
}
