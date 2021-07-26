#' plot_composition_heatmap
#'
#' @description Produces a heatmap showing the marker percentages within each cluster
#' and the cluster sizes
#' @param composition - a dataframe output from cluster_percent_composition
#' @param pheno_to_exclude Vector of phenotype to exclude
#' @param log_values TRUE if the percentages should be logged (base 10)
#' @param type_of_aggregate Cluster or Community
#' @param column Column with cell types
#' @importFrom grDevices colorRampPalette
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap anno_barplot
#' @importFrom reshape2 dcast
#' @return A plot is returned
#' @examples
#' communities <- identify_cell_communities(SPIAT::formatted_image, radius=100)
#' communities_vis <- composition_of_clusters_and_communities(communities, type_of_aggregate = "Community", column="Phenotype)
#' plot_composition_heatmap(communities_vis, type_of_aggregate = "Community", column="Phenotype")
#' @export

plot_composition_heatmap <- function(composition, pheno_to_exclude = NULL, log_values = FALSE, type_of_aggregate, column) {

  if(type_of_aggregate == "Community"){
    cluster_size <- unique(data.frame(Community = composition$Community,
                                      Total_cells = composition$Total_number_of_cells))
    rownames(cluster_size) <- cluster_size$Community
    cluster_size$Community <- NULL

    composition2 <- composition[,c(column, "Community", "Percentage")]
    composition2 <- dcast(composition2, paste(column, "~", type_of_aggregate), value.var="Percentage")

  }else if(type_of_aggregate == "Cluster"){
    cluster_size <- unique(data.frame(Cluster = composition$Cluster,
                                      Total_cells = composition$Total_number_of_cells))
    rownames(cluster_size) <- cluster_size$Cluster
    cluster_size$Cluster <- NULL

    composition2 <- composition[,c(column, "Cluster", "Percentage")]
    composition2 <- dcast(composition2, paste(column, "~", type_of_aggregate), value.var="Percentage")
  }else{
    stop("Only Community and Cluster are accepted as valid column names")
  }


  rownames(composition2) <- composition2[,column]
  composition2[,column] <- NULL
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
  map_cols <- colorRampPalette(c("white", "red"))(1000)
  ha = HeatmapAnnotation("size" = anno_barplot(cluster_size), show_annotation_name = FALSE)
  Heatmap(as.matrix(composition2), name=" ",
                          cluster_columns=TRUE,
                          cluster_rows=TRUE,
                          col = map_cols,
                          top_annotation = ha,
                          border = TRUE)
 
}