#' composition_of_clusters_and_communities
#'
#' @description Returns a dataframe which contains the percentages of cells with a
#' specific marker within each cluster and the number of cells in the cluster.
#' @param formatted_data_with_clusters - a dataframe output from generate_clusters
#' @param type_of_aggregate Cluster or Community
#' @param column Column with cell types
#' @importFrom stats aggregate
#' @return A data.frame is returned
#' @examples
#' communities <- identify_cell_communities(SPIAT::formatted_image, radius=100)
#' communities_vis <- composition_of_clusters_and_communities(communities, type_of_aggregate = "Community", column="Phenotype")
#' @export

composition_of_clusters_and_communities <- function(formatted_data_with_clusters, type_of_aggregate, column) {
    number_of_clusters <- length(unique(formatted_data_with_clusters[,type_of_aggregate]))

    colnames(formatted_data_with_clusters)[colnames(formatted_data_with_clusters) == column] <- "Temp_pheno"
    
    if(type_of_aggregate == "Community"){
      composition <- aggregate(Cell.ID ~ Temp_pheno + Community, formatted_data_with_clusters, length)
      colnames(composition)[3] <- "Number_of_cells"
      cluster_size <- table(formatted_data_with_clusters$Community)
      composition$Total_number_of_cells <- as.vector(cluster_size[match(composition$Community, names(cluster_size))])

    }else if(type_of_aggregate == "Cluster"){
      composition <- aggregate(Cell.ID ~ Temp_pheno + Cluster, formatted_data_with_clusters, length)
      colnames(composition)[3] <- "Number_of_cells"
      cluster_size <- table(formatted_data_with_clusters$Cluster)
      composition$Total_number_of_cells <- as.vector(cluster_size[match(composition$Cluster, names(cluster_size))])
    }else{
      stop("Only Community and Cluster are accepted as valid column names")
    }

    composition$Percentage <- (composition$Number_of_cells/composition$Total_number_of_cells)*100
    colnames(composition)[colnames(composition) == "Temp_pheno"] <- column
    return(composition)
}
