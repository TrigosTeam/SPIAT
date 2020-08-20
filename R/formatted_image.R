#' Example Inform image
#'
#' A dataset containing an example Inform image, loaded and formatted as described in data-raw/formatted_image.R 
#'
#' @format A SingleCellExperiment object where the count assay stores the expression level of every marker (rows) for
#every cell (columns), and cell phenotype, x and y coordinates, other properties (Cell Size, Nucleus Size, 
#' Nucleus Compactness, Nucleus Axis Ratio, Cell Axis Ratio ) are stored under colData
#under colData
#' @usage data(formatted_image)
"formatted_image"