#' Calculate the (normalised) mixing score for interested cell types
#'
#' @description Produces a data.frame with mixing scores of input reference and
#'   target cells from a SpatialExperiment object. It calculates
#'   reference-target interactions and reference-reference interactions based on
#'   a radius. It derives the mixing score and the normalised mixing score.
#'   Function returns NA if the mixing score is being calculated between cells
#'   of the same type.

#' @details The mixing score was originally defined as the number of
#'   immune-tumour interactions divided by the number of immune-immune
#'   interactions within a defined radius (Keren et al., 2018). The normalised
#'   mixing score normalises the immune-tumour interactions and immune-immune
#'   interactions within radius by the total number of immune-tumour and
#'   immune-immune interactions in the image, respectively. We have generalized
#'   this score to allow calculation of any two cell phenotypes defined by the
#'   user.
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param reference_celltype String Vector. Cell types of the reference cells.
#' @param target_celltype String Vector. Cell types of the target cells.
#' @param feature_colname String specifying the column with the desired cell
#'   type annotations.
#' @param radius The maximum radius around a reference cell type for another
#'   cell to be considered an interaction.
#' @import dplyr
#' @return A data.frame of cell numbers, number of cell interactions,  mixing
#'   scores, and normalised mixing scores. If there are no reference or target
#'   cells found in the image, or there are no reference cells found within the
#'   specified radius of any reference cells,the returned (normalised) mixing
#'   scores will be NA. If there are no target cells found within the radius of
#'   any refernece cells, the returned (normalised) mixing scores will be 0.
#' @export
#' @examples
#' mixing_score_summary(SPIAT::defined_image, reference_celltype = "Tumour", target_celltype="Immune1",
#' radius = 50, feature_colname = "Cell.Type")

mixing_score_summary <- function(spe_object, reference_celltype, target_celltype, 
                                 radius=20, feature_colname)
{
    formatted_data <- get_colData(spe_object)
    df.cols <- c("Reference", "Target", "Number_of_reference_cells",
                 "Number_of_target_cells", "Reference_target_interaction",
                 "Reference_reference_interaction", "Mixing_score", 
                 "Normalised_mixing_score")
    df <- data.frame(matrix(ncol=8,nrow=1, dimnames=list(NULL, df.cols)), stringsAsFactors = FALSE)
    for (i in reference_celltype) {
        reference_cells <- formatted_data[formatted_data[,feature_colname] == i,]
        for (j in target_celltype) {
            if (i == j) {
                df <-  rbind(df[ ,df.cols], 
                             c(i, j, nrow(reference_cells), nrow(target_cells), NA, NA, NA, NA))
            }   
            target_cells <- formatted_data[formatted_data[,feature_colname] == j,]
            if (nrow(reference_cells) == 0) {
                methods::show(paste("There are no unique reference cells of specified celltype", i, "for target cell", j))
                df <-  rbind(df[ ,df.cols], 
                             c(i, j, 0, nrow(target_cells), 0, 0, NA, NA))
            }
            else if (nrow(target_cells) == 0) {
                methods::show(paste("There are no unique target cells of specified celltype", j, "for reference cell", i))
                reference_cell_cords <- reference_cells[, c("Cell.X.Position", 
                                                            "Cell.Y.Position")]
                reference_reference_result <- dbscan::frNN(reference_cell_cords, 
                                                           eps = radius, sort = FALSE)
                reference_reference_interactions <- sum(rapply(reference_reference_result$id, 
                                                               length))/2 # halve it to avoid counting each ref-ref interaction twice
                df <-  rbind(df[ ,df.cols], 
                             c(i, j, nrow(reference_cells), 0, 0, reference_reference_interactions, NA, NA))
            }
            else{
                reference_cell_cords <- reference_cells[, c("Cell.X.Position", 
                                                            "Cell.Y.Position")]
                target_cell_cords <- target_cells[, c("Cell.X.Position", 
                                                      "Cell.Y.Position")]
                reference_target_result <- dbscan::frNN(target_cell_cords, eps = radius, 
                                                        query = reference_cell_cords, sort = FALSE)
                reference_target_interactions <- sum(rapply(reference_target_result$id, 
                                                            length))
                reference_reference_result <- dbscan::frNN(reference_cell_cords, 
                                                           eps = radius, sort = FALSE)
                reference_reference_interactions <- sum(rapply(reference_reference_result$id, 
                                                               length))/2 # halve it to avoid counting each ref-ref interaction twice
                if (reference_reference_interactions != 0) {
                    mixing_score <- reference_target_interactions/reference_reference_interactions
                    normalised_mixing_score <- 0.5 * mixing_score * (nrow(reference_cells)-1) / nrow(target_cells)
                }else {
                    normalised_mixing_score <- mixing_score <- 0
                    methods::show(paste("There are no reference to reference interactions for", j, "in the specified radius, cannot calculate mixing score"))
                    #stop()
                }
                df <-  rbind(df[ ,df.cols], 
                             c(i, j, nrow(reference_cells), nrow(target_cells), reference_target_interactions, reference_reference_interactions, mixing_score, normalised_mixing_score))
            }
        }
    }
    df[,3:8] <- vapply(df[,3:8],as.numeric,numeric(nrow(df[,3:8])))
    
    df <- df[-1,]
    return(df)
}
