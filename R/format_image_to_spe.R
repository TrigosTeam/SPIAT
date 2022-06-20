#' Format an image into a SpatialExperiment object
#'
#' @description Reads in spatial data in the form of cell coordinates, cell
#'   phenotypes (if available), and marker intensities and transforms to a
#'   SpatialExperiment object. The assay stores the intensity level of every
#'   marker (rows) for every cell (columns). Cell phenotype is stored under
#'   `colData()`. Cell x and y coordinates are stored under `spatialCoords()`
#'   field. 
#'
#' @export
#' @param intensity_matrix A matrix of marker intensities or gene expression
#'   where the column names are the Cell IDs, and the rownames the marker.
#' @param phenotypes (Optional) String Vector of cell phenotypes in the same
#'   order in which they appear in `intensity_matrix`. If no phenotypes
#'   available, then a vector of NAs can be used as input. Note that the
#'   combination of markers (e.g. CD3,CD4) needs to be used instead of the cell
#'   type name (e.g. helper T cells).
#' @param coord_x Numeric Vector with the X coordinates of the cells. The cells
#'   must be in the same order as in the `intensity_matrix`.
#' @param coord_y Numeric Vector with the Y coordinates of the cells. The cells
#'   must be in the same order as in the `intensity_matrix`.
#' @return A SpatialExperiment object is returned
#' @examples
#' #Construct a marker intensity matrix (rows are markers, columns are cells)
#' intensity_matrix <- matrix(c(14.557, 0.169, 1.655, 0.054, 17.588, 0.229,
#' 1.188, 2.074, 21.262, 4.206,  5.924, 0.021), nrow = 4, ncol = 3)
#' # define marker names as rownames
#' rownames(intensity_matrix) <- c("DAPI", "CD3", "CD4", "AMACR")
#' # define cell IDs as colnames
#' colnames(intensity_matrix) <- c("Cell_1", "Cell_2", "Cell_3")
#' # Construct a dummy metadata (phenotypes, x/y coordinates)
#' # the order of the elements in these vectors correspond to the cell order
#' # in `intensity matrix`
#' phenotypes <- c("OTHER",  "AMACR", "CD3,CD4")
#' coord_x <- c(82, 171, 184)
#' coord_y <- c(30, 22, 38)
#'
#' formatted_image <- format_image_to_spe(intensity_matrix=intensity_matrix,
#' phenotypes = phenotypes, coord_x = coord_x,coord_y = coord_y)

format_image_to_spe <- function(intensity_matrix, phenotypes = NULL, coord_x, 
                                coord_y){
  
  intensity_matrix <- remove_intensity_na(intensity_matrix)

  metadata_columns <- data.frame(Phenotype = phenotypes,
                                 Cell.X.Position = coord_x,
                                 Cell.Y.Position = coord_y)
  spe <- SpatialExperiment::SpatialExperiment(
    assay = intensity_matrix,
    colData = metadata_columns,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position"))
    
  return(spe)
}
