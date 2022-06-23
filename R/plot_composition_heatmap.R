#' plot_composition_heatmap
#'
#' @description Produces a heatmap showing the marker percentages within each
#'   cluster and the cluster sizes.
#' @param composition Data.frame. Output from
#'   \code{\link{composition_of_neighborhoods}}.
#' @param pheno_to_exclude String Vector of phenotype to exclude.
#' @param log_values Boolean. TRUE if the percentages should be logged (base
#'   10).
#' @param feature_colname String. Column with cell types.
#' @return A plot is returned
#' @examples
#' neighborhoods <- identify_neighborhoods(image_no_markers, method = "hierarchical",
#' min_neighborhood_size = 100, cell_types_of_interest = c("Immune", "Immune1", "Immune2"),
#' radius = 50, feature_colname = "Cell.Type")
#' neighborhoods_vis <- composition_of_neighborhoods(neighborhoods, feature_colname="Cell.Type")
#' plot_composition_heatmap(neighborhoods_vis, feature_colname="Cell.Type")
#' @export

plot_composition_heatmap <- function(composition, pheno_to_exclude = NULL, 
                                     log_values = FALSE, feature_colname) {
    cluster_size <- unique(data.frame(
        Neighborhood = composition$Neighborhood,
        Total_cells = composition$Total_number_of_cells))
    rownames(cluster_size) <- cluster_size$Neighborhood
    cluster_size$Neighborhood <- NULL
    
    composition2 <- composition[,c(feature_colname, 
                                   "Neighborhood", "Percentage")]
    composition2 <- reshape2::dcast(composition2, 
                                    paste(feature_colname, "~", "Neighborhood"),
                                    value.var="Percentage")
    
    rownames(composition2) <- composition2[,feature_colname]
    composition2[,feature_colname] <- NULL
    composition2[is.na(composition2)] <- -1
    composition2 <- as.matrix(composition2)
    
    if(!is.null(pheno_to_exclude)){
        composition2 <- composition2[!(rownames(composition2) %in% 
                                           pheno_to_exclude),]
    }
    
    if(log_values){
        composition2 <- apply(composition2, 2, log10)
        composition2[is.nan(composition2)] <- 1000
        composition2[composition2 == 1000] <- min(composition2)-1
    }
    
    #plot the heatmap
    map_cols <- grDevices::colorRampPalette(c("white", "red"))(1000)
    ha <- ComplexHeatmap::HeatmapAnnotation(
        "size" = ComplexHeatmap::anno_barplot(cluster_size), 
        show_annotation_name = FALSE)
    ComplexHeatmap::Heatmap(as.matrix(composition2), name=" ",
                            cluster_columns=TRUE,
                            cluster_rows=TRUE,
                            col = map_cols,
                            top_annotation = ha,
                            border = TRUE)
    
}
