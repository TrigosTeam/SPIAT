#' select_celltypes
#'
#' @description Select cell types to keep or exclude in the analysis. The output
#'   of this function also includes the original image size and cell count.
#'
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param celltypes String Vector of celltypes of keep or exclude.
#' @param feature_colname String. The column that has the interested cell types.
#'   If the cells ids are used to select cells, use "Cell.ID" for this arg.
#' @param keep Boolean. TRUE if vector of `celltypes` are the cells that are
#'   going to be kept, FALSE if they are to be removed.
#' @return A SpatialExperiment object is returned. The original image size
#'   and cell count can be accessed by `attr(slim_spe, "original_cell_number")`
#'   and `attr(slim_spe, "range_of_coords")`, where `slim_spe` is the output of
#'   this function.
#'
#' @examples
#' data_subset <- select_celltypes(SPIAT::simulated_image,
#' celltypes = c("Tumour_marker","Immune_marker1","Immune_marker2",
#' "Immune_marker3","Immune_marker4"),
#' feature_colname = "Phenotype", keep=TRUE)
#' attr(data_subset, "original_cell_number") #cell number in the original image
#' attr(data_subset, "range_of_coords")
#' dim(data_subset)[2] # this is the new image cell number
#' @export
select_celltypes <- function(spe_object, celltypes, 
                             feature_colname = "Phenotype", keep = TRUE){
    data <- get_colData(spe_object)
    
    # # if rownames have the cell types of interest, make rownames as a column
    if (feature_colname == "Cell.ID"){
        SummarizedExperiment::colData(spe_object)$Cell.ID <- data$Cell.ID
    }
    
    # remember the total number of cells
    n_cells <- dim(data)[1]
    # remember the range of the cell coordinates
    xmax <- max(data$Cell.X.Position)
    xmin <- min(data$Cell.X.Position)
    ymax <- max(data$Cell.Y.Position)
    ymin <- min(data$Cell.Y.Position)
    
    # delete the cell rows
    if (keep) slim_spe <- 
        spe_object[, (spe_object[[feature_colname]] %in% celltypes)]
    else slim_spe <- 
        spe_object[, !(spe_object[[feature_colname]] %in% celltypes)]
    
    # # delete the rownames from the spe_object 
    # SummarizedExperiment::colData(spe_object)$rownames <- NULL
    
    # save the original image info in the slim attr
    attr(slim_spe, "original_cell_number") <- n_cells
    attr(slim_spe, "range_of_coords") <- c(xmin = xmin, xmax = xmax, 
                                           ymin = ymin, ymax = ymax)
    
    return(slim_spe)
}
