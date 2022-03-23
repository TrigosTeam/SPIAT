#' mixing_score_summary
#'
#' @description Produces a data.frame with mixing scores of inputed reference and target cells from a SingleCellExperiment object. 
#' It calculates reference-target interactions and reference-reference interactions based on default radius of 20.
#' It derives the mixing score and normalises the score by (mixing score) * (number of reference cells) * 2 / (number of target cells).
#' Function returns NA if the mixing score is being calculated between cells of the same type --include reference for mixing score itself?
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param reference_celltype Cell types of the reference cells
#' @param target_celltype Cell types of the target cells
#' @param feature_colname String specifying the column with the desired cell type annotations 
#' @param radius The maximum radius around a reference celltype for another cell to be considered an interaction.
#' @import dplyr
#' @return A data.frame of cell numbers, mixing scores, and normalised mixing scores.
#' @export
#' @examples 
#' mixing_score_summary(SPIAT::defined_image, reference_celltype = "Tumour", target_celltype="Immune1",
#' radius = 50, feature_colname = "Cell.Type")

mixing_score_summary <- function(sce_object, reference_celltype, target_celltype, 
                                 radius=20, feature_colname)
{
    formatted_data <- get_colData(sce_object)
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
            tryCatch({
                target_cells <- formatted_data[formatted_data[,feature_colname] == j,]
                if (nrow(reference_cells) == 0) {
                    print(paste("There are no unique reference cells of specified celltype", i, "for target cell", j))
                }
                if (nrow(target_cells) == 0) {
                    print(paste("There are no unique target cells of specified celltype", j, "for reference cell", i))
                }
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
                                                               length))
                if (reference_reference_interactions != 0) {
                    mixing_score <- reference_target_interactions/reference_reference_interactions
                    normalised_mixing_score <- 2 * mixing_score * (nrow(reference_cells)-1) / nrow(target_cells)
                }else {
                    normalised_mixing_score <- mixing_score <- NA
                    print(paste("There are no reference to reference interactions for", j, "in the specified radius, cannot calculate mixing score"))
                    #stop()
                }
                df <-  rbind(df[ ,df.cols], 
                             c(i, j, nrow(reference_cells), nrow(target_cells), reference_target_interactions, reference_reference_interactions, mixing_score, normalised_mixing_score))
            }, error=function(e){})
        }
    }
    df[,3:8] <- sapply(df[,3:8],as.numeric)
    df <- df[-1,]
    return(df)
}
