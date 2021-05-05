#' identify_cell_communities
#'
#' @description Identifies communities of cells based on their location. It excludes cells without a phenotype
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param clustering_method String specifying which clustering algorithm to use. Current options:
#' "dbscan" and "rphenograph"
#' @param radius Integer specifying the radius of search. Required for "dbscan"
#' @param min_community_size Minimum number of cells in a community
#' @param phenotypes_of_interest Vector of phenotypes to consider
#' @import ggplot2
#' @import dplyr
#' @import Rphenograph
#' @importFrom SummarizedExperiment colData assay
#' @importFrom tibble rownames_to_column
#' @importFrom dbscan dbscan
#' @importFrom dittoSeq dittoColors
#' @return A data.frame and a plot is returned
#' @examples
#' @export

identify_cell_communities <- function(sce_object, clustering_method = "dbscan", radius = NULL, min_community_size = 50, phenotypes_of_interest = NULL){

    # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
    Cell.X.Position <- Cell.Y.Position <- Community <- Xpos <- Ypos <- community <- NULL
      
    formatted_data <- data.frame(colData(sce_object))
    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

    intensity_matrix <- assay(sce_object)

    markers <- rownames(intensity_matrix)
    cell_ids <- colnames(intensity_matrix)

    rownames(intensity_matrix) <- NULL
    colnames(intensity_matrix) <- NULL
    intensity_matrix_t <- t(intensity_matrix)
    intensity_df <- data.frame(intensity_matrix_t)
    colnames(intensity_df) <- markers

    formatted_data <- cbind(formatted_data, intensity_df)
    formatted_data <- formatted_data[complete.cases(formatted_data),]

    ######remove cells without a phenotype
    formatted_data <- formatted_data[formatted_data$Phenotype != "OTHER", ]
    cell_cords <- formatted_data[,c("Cell.X.Position", "Cell.Y.Position")]

    #select cells to include if phenotypes of interest are specified
    if (!is.null(phenotypes_of_interest)) {
        formatted_data <- formatted_data[formatted_data$Phenotype %in% phenotypes_of_interest,]
    }
    
    #CHECK
    if (nrow(formatted_data) == 0) {
      stop("There are no cells in data/no cells for the phenotypes of interest")
    }

    if (clustering_method == "rphenograph") {
      Rphenograph_out <- Rphenograph(cell_cords, k = min_community_size)
      formatted_data$Community <- factor(membership(Rphenograph_out[[2]]))
    } else if (clustering_method == "dbscan") {
        #Use dbscan to generate clusters
        db <- dbscan::dbscan(cell_cords, eps = radius, minPts = min_community_size)
        #since dbscan outputs cluster 0 as noise, we add 1 to all cluster numbers to keep it consistent
        formatted_data$Community <- factor(db$cluster + 1)
    } else {
        stop("Please select a valid clustering method, current options: dbscan")
    }

    #start a plot for visualizing communities
    q <- ggplot(formatted_data, aes(x=Cell.X.Position, y=Cell.Y.Position))
    q <- q + geom_point(aes(color = Community), size = 0.01)

    #get number_of_communities
    number_of_communities <- length(unique(formatted_data$Community))
    
    # use colourblind-friendly colours
    colours <- dittoColors()[seq_len(number_of_communities)]

    #label the community centre by averaging x and y
    label_location <- vector()
    for (communitynumber in seq_len(number_of_communities)) {
        cells_in_community <- formatted_data[formatted_data$Community == communitynumber, ]
        minX <- min(cells_in_community$Cell.X.Position)
        maxX <- max(cells_in_community$Cell.X.Position)
        minY <- min(cells_in_community$Cell.Y.Position)
        maxY <- max(cells_in_community$Cell.Y.Position)
        averageX <- (minX + maxX)/2
        averageY <- (minY + maxY)/2

        label_location <- rbind(label_location,c(communitynumber, averageX, averageY))
    }
    label_location <- as.data.frame(label_location)
    colnames(label_location) <- c("community", "Xpos", "Ypos")
    q <- q + geom_text(data = label_location, aes(x = Xpos, y = Ypos, label = community))

    #label x and y axis and community labels
    q <- q + xlab("Cell.X.Position") + ylab("Cell.Y.Position")
    q <- q + scale_color_manual(values=colours)
    q <- q + theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.position = "none")
    print(q)

    formatted_data_with_communities <- formatted_data
    formatted_data_with_communities$Community <- paste0("Community_", formatted_data_with_communities$Community)
    return(formatted_data_with_communities)
}

