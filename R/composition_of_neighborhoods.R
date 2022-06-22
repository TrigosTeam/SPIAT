#' composition_of_neighborhoods
#'
#' @description Returns a data.frame which contains the percentages of cells
#'   with a specific marker within each neighborhood. and the number of cells in
#'   the neighborhood.
#' @param spe_object SpatialExperiment that is the output of
#'   \code{\link{identify_neighborhoods}}.
#' @param feature_colname String. Column with cell types.
#' @return A data.frame is returned
#' @examples
#' neighborhoods <- identify_neighborhoods(image_no_markers,
#' method = "hierarchical", min_neighborhood_size = 100,
#' cell_types_of_interest = c("Immune", "Immune1", "Immune2"), radius = 50,
#' feature_colname = "Cell.Type")
#' neighborhoods_vis <- composition_of_neighborhoods(neighborhoods,
#' feature_colname="Cell.Type")
#' @export

composition_of_neighborhoods <- function(spe_object, feature_colname) {
    
    neighborhoods_df <- get_colData(spe_object)
    neighborhoods_df <- 
        neighborhoods_df[stats::complete.cases(neighborhoods_df),]
    
    number_of_clusters <- length(unique(neighborhoods_df[,"Neighborhood"]))

    colnames(neighborhoods_df)[colnames(neighborhoods_df) == feature_colname] <-
        "Temp_pheno"
    
    composition <- stats::aggregate(Cell.ID ~ Temp_pheno + Neighborhood, 
                                    neighborhoods_df, length)
    colnames(composition)[3] <- "Number_of_cells"
    cluster_size <- table(neighborhoods_df$Neighborhood)
    composition$Total_number_of_cells <- 
        as.vector(cluster_size[match(composition$Neighborhood, 
                                     names(cluster_size))])

    composition$Percentage <- 
        (composition$Number_of_cells/composition$Total_number_of_cells)*100
    colnames(composition)[colnames(composition) == "Temp_pheno"] <- 
        feature_colname
    return(composition)
}
