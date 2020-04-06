#' composition_of_clusters_and_communities
#'
#' @description Returns a dataframe which contains the percentages of cells with a
#' specific marker within each cluster and the number of cells in the cluster.
#' @param formatted_data_with_clusters - a dataframe output from generate_clusters
#' @param column_to_consider Column name to consider as community/clusters
#' @export

composition_of_clusters_and_communities <- function(formatted_data_with_clusters, column_to_consider) {
    number_of_clusters <- length(unique(formatted_data_with_clusters[,column_to_consider]))

    if(column_to_consider == "Community"){
      composition <- aggregate(Cell.ID ~ Phenotype + Community, formatted_data_with_clusters, length)
      colnames(composition)[3] <- "Number_of_cells"
      cluster_size <- table(formatted_data_with_clusters$Community)
      composition$Total_number_of_cells <- as.vector(cluster_size[match(composition$Community, names(cluster_size))])

    }else if(column_to_consider == "Cluster"){
      composition <- aggregate(Cell.ID ~ Phenotype + Cluster, formatted_data_with_clusters, length)
      colnames(composition)[3] <- "Number_of_cells"
      cluster_size <- table(formatted_data_with_clusters$Cluster)
      composition$Total_number_of_cells <- as.vector(cluster_size[match(composition$Cluster, names(cluster_size))])
    }else{
      print("Only Community and Cluster are accepted as valid column names")
    }

    composition$Percentage <- (composition$Number_of_cells/composition$Total_number_of_cells)*100
    return(composition)
}
