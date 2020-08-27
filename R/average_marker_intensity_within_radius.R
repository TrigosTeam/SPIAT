#' average_marker_intensity_within_radius
#'
#' @description Calculates the average intensity of the target_marker within a radius
#' from the cells positive for the reference marker.
#' Note that it pools all cells with the target marker that are within the specific radius of
#' any reference cell. Results represent the average intensities within a radius,
#' but do not correspond to metrics for each cell
#'
#' @param sce_object A SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param reference_marker A string specifying the marker that is used for reference cells
#' @param target_marker A string specifying the marker to calculate its average intensity
#' @param radius An integer specifying the radius of search for cells around the reference cells
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom dbscan frNN
#' @importFrom stats complete.cases
#' @importFrom SummarizedExperiment assay colData
#' @export

average_marker_intensity_within_radius <- function(sce_object, reference_marker, target_marker, radius = 20) {

    formatted_data <- data.frame(colData(sce_object))
    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

    intensity_matrix <- assay(sce_object)

    markers <- rownames(intensity_matrix)
    cell_ids <- colnames(intensity_matrix)

    rownames(intensity_matrix) <- NULL
    colnames(intensity_matrix) <- NULL
    intensity_matrix_t <- t(intensity_matrix)
    intensity_df <- data.frame(intensity_matrix_t)
    colnames(intensity_df) <- markers

    formatted_data <- cbind(formatted_data, intensity_df)
    formatted_data <- formatted_data[complete.cases(formatted_data), ]

    #Select the cells that have the reference marker phenotype
    reference_cells <- formatted_data[grepl(reference_marker, formatted_data$Phenotype),]
    if (nrow(reference_cells) == 0) {
        stop("There are no reference cells found for the marker")
    }
    
    #Target cells are don't contain the reference marker
    target_cells <- formatted_data[grepl(target_marker, formatted_data$Phenotype),]
    if (nrow(target_cells) == 0) {
        stop("There are no target cells found for the marker")
    }
    
    #Remove cells with both markers
    common_cells <- reference_cells$Cell.ID[reference_cells$Cell.ID %in% target_cells$Cell.ID]

    reference_cells <- reference_cells[!(reference_cells$Cell.ID %in% common_cells),]
    target_cells <- target_cells[!(target_cells$Cell.ID %in% common_cells),]

    #Get the coordinates to find neighbours
    reference_cell_cords <- reference_cells[,c("Cell.X.Position", "Cell.Y.Position")]
    target_cell_cords <- target_cells[,c("Cell.X.Position", "Cell.Y.Position")]

    #frNN output ids, the rowid of reference_cell_cords matches the row number of target_cell_cords
    search_result <- frNN(target_cell_cords, eps = radius, query = reference_cell_cords, sort = FALSE)
    rownums <- unique(unlist(search_result$id))

    #check
    if (length(rownums) == 0) {
        stop("There are no target cells within the specified radius, cannot calculate average intensity")
    } else {
        target_within_radius <- target_cells[rownums,]
        average_marker_intensity <- mean(target_within_radius[,target_marker], na.rm=TRUE)
    }
    return(average_marker_intensity)
}
