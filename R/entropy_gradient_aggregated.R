#' The aggregated gradient of entropy and the peak of the gradient
#'
#' @description This function first calculates the entropy within circles of
#'   each reference cell at each radius. Then at each radius, the entropy of all
#'   circles surrounding each cell are aggregated into one number. The function
#'   sweeps over the specified radii and calculates the aggregated entropy under
#'   each radius.
#'
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param cell_types_of_interest String Vector. The cell types that the entropy
#'   is computed on. 
#' @param feature_colname String. The column name of the interested cell types.
#' @param radii Numeric Vector. A vector of radii within a circle of a reference
#'   cell where the entropy is computed on.
#' @import dplyr
#'
#' @return A list of the gradient of entropy and the peak
#' @export
#'
#' @examples
#' gradient_pos <- seq(50, 500, 50)
#' gradient_results <- entropy_gradient_aggregated(SPIAT::defined_image,
#' cell_types_of_interest = c("Tumour","Immune3"),
#' feature_colname = "Cell.Type", radii = gradient_pos)
#' plot(1:10,gradient_results$gradient_df[1, 3:12])
entropy_gradient_aggregated <- function(spe_object, cell_types_of_interest,
                                        feature_colname, radii){
    
    # entropy_pairs <- vector()
    gradient_aggregated <- vector()
    
    for(pheno1 in cell_types_of_interest){
        for(pheno2 in cell_types_of_interest){
            if(pheno2 != pheno1){
                unique_types <- unique(SummarizedExperiment::colData(spe_object)
                                       [,feature_colname])
                if((pheno1 %in% unique_types) && (pheno2 %in% unique_types)){
                    # get the gradient
                    gradient_local <- compute_gradient(
                        spe_object, radii = radii,number_of_cells_within_radius,
                        reference_celltype = pheno1,
                        target_celltype = c(pheno1, pheno2),
                        feature_colname = feature_colname)
                    
                    
                    gradient_local_aggregated <- vector()
                    # aggregate all tumour and immune counts and calculate new aggragated "entropy"
                    for(j in seq_len(length(gradient_local))){
                        temp <- gradient_local[[j]]
                        total_reference <- sum(temp[[pheno1]][[pheno1]])
                        total_target <- sum(temp[[pheno1]][[pheno2]])
                        total_cells <- total_reference + total_target
                        reference_disorder <- (total_reference/total_cells) * 
                            log2(total_reference/total_cells)
                        target_disorder <- (total_target/total_cells) * 
                            log2(total_target/total_cells)
                        local_entropy <- -(reference_disorder+target_disorder)
                        gradient_local_aggregated <- 
                            c(gradient_local_aggregated,local_entropy)
                        
                    }
                    gradient_aggregated <- rbind(
                        gradient_aggregated,c(pheno1, pheno2, 
                                              gradient_local_aggregated))
                }
                else{
                    gradient_local_aggregated <- 
                        c(pheno1, pheno2, rep(NA, length(radii)))
                    gradient_aggregated <- rbind(gradient_aggregated,
                                                 gradient_local_aggregated)
                }}}}
    
    # format the output
    colnames(gradient_aggregated) <- c("Celltype1", "Celltype2",
                                       paste("Pos", radii, sep="_"))
    gradient_aggregated <- as.data.frame(gradient_aggregated)
    gradient_aggregated[,3:(length(radii) + 2)] <-
        apply(gradient_aggregated[,3:(length(radii) + 2)], 2,
              function(x){as.numeric(as.character(x))})
    rownames(gradient_aggregated) <- NULL
    
    # get the peak of the gradient
    temp <- gradient_aggregated[1,]
    temp$Celltype1 <- NULL
    temp$Celltype2 <- NULL
    temp <- as.numeric(temp)
    peak <- which(temp == max(temp, na.rm = TRUE))
    
    return(list(gradient_df = gradient_aggregated, peak = peak))
}
