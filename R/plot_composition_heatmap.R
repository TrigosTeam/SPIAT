#' plot_composition_heatmap
#'
#' @description Produces a heatmap showing the marker percentages within each cluster
#' and the cluster sizes
#' @param composition - a dataframe output from cluster_percent_composition
#' @param pheno_to_exclude Vector of phenotype to exclude
#' @param log_values TRUE if the percentages should be logged (base 10)
#' @param column_to_consider Column name to consider as community/clusters
#' @import RColorBrewer
#' @importFrom pheatmap pheatmap
#' @importFrom reshape2 dcast
#' @export

#ColorRampPalette is from 'dichromat' package but loaded by 'pheatmap'

plot_composition_heatmap <- function(composition, pheno_to_exclude = NULL, log_values = FALSE, column_to_consider) {

  if(column_to_consider == "Community"){
    cluster_size <- unique(data.frame(Community = composition$Community,
                                      Total_cells = composition$Total_number_of_cells))
    rownames(cluster_size) <- cluster_size$Community
    cluster_size$Community <- NULL

    composition2 <- composition[,c("Phenotype", "Community", "Percentage")]
    composition2 <- dcast(composition2,  Phenotype ~ Community, value.var="Percentage")

  }else if(column_to_consider == "Cluster"){
    cluster_size <- unique(data.frame(Cluster = composition$Cluster,
                                      Total_cells = composition$Total_number_of_cells))
    rownames(cluster_size) <- cluster_size$Cluster
    cluster_size$Cluster <- NULL

    composition2 <- composition[,c("Phenotype", "Cluster", "Percentage")]
    composition2 <- dcast(composition2,  Phenotype ~ Cluster, value.var="Percentage")
  }else{
    stop("Only Community and Cluster are accepted as valid column names")
  }


  rownames(composition2) <- composition2$Phenotype
  composition2$Phenotype <- NULL
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
    composition2[is.nan(composition2_log10)] <- 1000
    composition2[composition2_log10 == 1000] <- min(composition2_log10)-1
  }

  #plot the heatmap
  map_cols <- colorRampPalette(c("white", "red"))(100)
  anno_cols <- list(Total_cells = c("white", "blue"))
  pheatmap(as.matrix(composition2), color = map_cols,
           #annotation_col = cluster_size, 
           annotation_colors = anno_cols,
           cluster_cols=TRUE, cluster_rows=TRUE)
}
