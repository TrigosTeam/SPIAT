#' composition_of_neighborhoods
#'
#' @description Returns a dataframe which contains the percentages of cells with
#'   a specific marker within each cluster and the number of cells in the
#'   cluster.
#' @param formatted_data_with_clusters - a dataframe output from
#'   generate_clusters
#' @param feature_colname Column with cell types
#' @return A data.frame is returned
#' @examples
#' communities <- identify_neighborhoods(SPIAT::formatted_image, radius=100)
#' communities_vis <- composition_of_neighborhoods(communities, feature_colname="Phenotype")
#' @export

composition_of_neighborhoods <- function(formatted_data_with_clusters, feature_colname) {
    number_of_clusters <- length(unique(formatted_data_with_clusters[,"Cluster"]))

    colnames(formatted_data_with_clusters)[colnames(formatted_data_with_clusters) == feature_colname] <- "Temp_pheno"
    
    composition <- stats::aggregate(Cell.ID ~ Temp_pheno + Cluster, formatted_data_with_clusters, length)
    colnames(composition)[3] <- "Number_of_cells"
    cluster_size <- table(formatted_data_with_clusters$Cluster)
    composition$Total_number_of_cells <- as.vector(cluster_size[match(composition$Cluster, names(cluster_size))])

    composition$Percentage <- (composition$Number_of_cells/composition$Total_number_of_cells)*100
    colnames(composition)[colnames(composition) == "Temp_pheno"] <- feature_colname
    return(composition)
}
