#' calculate the distances of each cell to the tumour margin
#'
#' @description Returns a SPE object with the minimum distance from immune cells
#'   to the identified tumour bordering cells.
#'
#' @param spe_object SpatialExperiment object. It should contain information of
#'   the detected bordering cells (`colData()` has `Region` column).
#' @return An spe_object is returned
#' @export
#' @examples
#' spe_border <- identify_bordering_cells(SPIAT::defined_image, 
#' reference_cell = "Tumour", feature_colname = "Cell.Type", n_to_exclude = 10)
#' spe_dist <- calculate_distance_to_tumour_margin(spe_border)

calculate_distance_to_tumour_margin <- function(spe_object){
    
    #CHECK if the user has found the bordering cells yet
    if (is.null(spe_object$Region)){
        stop("Please find the bordering cells first! (use identify_bordering_cells)")
    }
    
    dat <- get_colData(spe_object)
    
    #CHECK
    if (nrow(dat) == 0) {
        stop("There are no cells")
    }
    
    dat <- dat[, c("Cell.ID","Region", "Cell.X.Position", "Cell.Y.Position")]
    dat <- dat[dat$Region != "",]
    
    spe_object <- define_celltypes(
        spe_object, categories = c("Outside","Inside","Border"),
        category_colname ="Region", 
        names = c("Non-border","Non-border","Border"), new_colname = "region2")
    
    dist_matrix <- calculate_minimum_distances_between_celltypes(
        spe_object, cell_types_of_interest = c("Non-border","Border"), 
        feature_colname ="region2")
    dist_matrix[dist_matrix$RefType == "Border", "Distance"] <- 0
    dist_matrix$Order <- as.numeric(substr(dist_matrix$RefCell,
                                           start = 6, stop = 30))
    dist_matrix <- dist_matrix[order(dist_matrix$Order),]
    dat <- merge(dat, dist_matrix,by.x = "Cell.ID",by.y = "RefCell", 
                 all.x = TRUE, sort = FALSE)
    SummarizedExperiment::colData(spe_object)$Distance.To.Border <- dat$Dist
    
    return(spe_object)
}
