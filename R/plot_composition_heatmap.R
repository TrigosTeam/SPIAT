#' plot_composition_heatmap
#'
#' @description Produces a heatmap showing the marker percentages within each cluster
#' and the cluster sizes
#' @param composition - a dataframe output from cluster_percent_composition
#' @param pheno_to_exclude Vector of phenotype to exclude
#' @param log_values TRUE if the percentages should be logged (base 10)
#' @param type_of_aggregate Cluster or Community
#' @param feature_colname String. Column with cell types.
#' @return A plot is returned
#' @examples
#' communities <- identify_cell_communities(SPIAT::formatted_image, radius=100)
#' communities_vis <- composition_of_clusters_and_communities(communities, type_of_aggregate = "Community", feature_colname="Phenotype)
#' plot_composition_heatmap(communities_vis, type_of_aggregate = "Community", feature_colname="Phenotype")
#' @export

plot_composition_heatmap <- function(composition, pheno_to_exclude = NULL, log_values = FALSE, type_of_aggregate, feature_colname) {

  if(type_of_aggregate == "Community"){
    cluster_size <- unique(data.frame(Community = composition$Community,
                                      Total_cells = composition$Total_number_of_cells))
    rownames(cluster_size) <- cluster_size$Community
    cluster_size$Community <- NULL

    composition2 <- composition[,c(feature_colname, "Community", "Percentage")]
    composition2 <- reshape2::dcast(composition2, paste(feature_colname, "~", type_of_aggregate), value.var="Percentage")

  }else if(type_of_aggregate == "Cluster"){
    cluster_size <- unique(data.frame(Cluster = composition$Cluster,
                                      Total_cells = composition$Total_number_of_cells))
    rownames(cluster_size) <- cluster_size$Cluster
    cluster_size$Cluster <- NULL

    composition2 <- composition[,c(feature_colname, "Cluster", "Percentage")]
    composition2 <- reshape2::dcast(composition2, paste(feature_colname, "~", type_of_aggregate), value.var="Percentage")
  }else{
    stop("Only Community and Cluster are accepted as valid colname names")
  }


  rownames(composition2) <- composition2[,feature_colname]
  composition2[,feature_colname] <- NULL
  composition2[is.na(composition2)] <- -1
  composition2 <- as.matrix(composition2)

  if(!is.null(pheno_to_exclude)){
    composition2 <- composition2[!(rownames(composition2) %in% pheno_to_exclude),]
    ##Remove communities that are only zeros
    #total_NA <- apply(composition2, 2, function(x){ sum(is.na(x))})
    #max_NA <- nrow(composition2)
    #com_to_remove <- names(total_NA)[total_NA == max_NA]
    #composition2 <- composition2[,!(colnames(composition2) %in% com_to_remove)]
  }

  if(log_values){
    composition2 <- apply(composition2, 2, log10)
    composition2[is.nan(composition2)] <- 1000
    composition2[composition2 == 1000] <- min(composition2)-1
  }

  #plot the heatmap
  map_cols <- grDevices::colorRampPalette(c("white", "red"))(1000)
  ha = ComplexHeatmap::HeatmapAnnotation("size" = ComplexHeatmap::anno_barplot(cluster_size), 
                                         show_annotation_name = FALSE)
  ComplexHeatmap::Heatmap(as.matrix(composition2), name=" ",
                          cluster_columns=TRUE,
                          cluster_rows=TRUE,
                          col = map_cols,
                          top_annotation = ha,
                          border = TRUE)
 
}