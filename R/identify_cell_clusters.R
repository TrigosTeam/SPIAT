#' identify_cell_clusters
#'
#' @description Uses Euclidean distances to identify clusters of cells within a specified radius.
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param cell_types_of_interest Vector of phenotypes to consider
#' @param radius Integer specifying the radius of search.
#' @param column Column from which the cell types are selected
#' @param no_pheno Cell type corresponding to cells without a known phenotype (e.g. "None", "Other")
#' @import dplyr
#' @importFrom SummarizedExperiment colData assay
#' @importFrom tibble rownames_to_column
#' @importFrom stats complete.cases hclust cutree as.dist
#' @importFrom apcluster negDistMat
#' @importFrom dittoSeq dittoColors
#' @import ggplot2
#' @return A data.frame and a plot is returned
#' @examples
#' @export

# imported ggplo2 as interdependency of functions

identify_cell_clusters <- function(sce_object, cell_types_of_interest, radius, column, no_pheno = NULL) {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Cell.X.Position <- Cell.Y.Position <- Cluster <- Xpos <- Ypos <- NULL
  
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
  if(!is.null(no_pheno)){
    formatted_data <- formatted_data[formatted_data[,column] != no_pheno, ]
  }

  #select cells to include if phenotypes of interest are specified
  if (!is.null(cell_types_of_interest)) {
    formatted_data <- formatted_data[formatted_data[,column] %in% cell_types_of_interest,]
  }
  
  #CHECK
  if (nrow(formatted_data) == 0) {
    stop("There are no cells in data/no cells for the phenotypes of interest")
  }

  rownames(formatted_data) <- formatted_data$Cell.ID
  sim_close <- - negDistMat(formatted_data[,c("Cell.X.Position", "Cell.Y.Position")])

  sim_close[sim_close > radius] <- NA
  sim_close[sim_close == 0] <- NA

  #Remove rows, columns with only NAs
  sim_close <- sim_close[apply(sim_close, 1, function(x){
      if(sum(is.na(x)) == length(x)){
        return(FALSE)
      }else{
        return(TRUE)
      }} ),]

  if(!is.null(dim(sim_close))){
      sim_close <- sim_close[,apply(sim_close, 2, function(x){
        if(sum(is.na(x)) == length(x)){
          return(FALSE)
        }else{
          return(TRUE)
        }} )]

      cells_in_cluster <- rownames(sim_close)
      sim_close <- ifelse(is.na(sim_close), 1, 0)
      if(nrow(sim_close) != 0 & ncol(sim_close) != 0){
        h <- hclust(as.dist(sim_close), method="single")

        local_clusters <- cutree(h, h = 0.5)

        formatted_data$Cluster <- as.character(local_clusters[match(formatted_data$Cell.ID, names(local_clusters))])
        
      } else {
        formatted_data$Cluster <- NA
        #stop("The radius specified may be too small, no clusters were found")
      }
  } else {
    formatted_data$Cluster <- NA  
    #stop("The radius specified may be too small, no clusters were found")
  }

  #get cells assigned to clusters
  cells_in_clusters <- formatted_data[complete.cases(formatted_data),]
  cells_not_in_clusters <- formatted_data[!complete.cases(formatted_data),]
  
  #get number_of_clusters
  number_of_clusters <- length(unique(cells_in_clusters$Cluster))

  #label the Cluster centre by averaging x and y

  label_location <- vector()
  for (Clusternumber in seq_len(number_of_clusters)) {
    cells_in_Cluster <- cells_in_clusters[cells_in_clusters$Cluster == Clusternumber, ]
    minX <- min(cells_in_Cluster$Cell.X.Position)
    maxX <- max(cells_in_Cluster$Cell.X.Position)
    minY <- min(cells_in_Cluster$Cell.Y.Position)
    maxY <- max(cells_in_Cluster$Cell.Y.Position)
    averageX <- (minX + maxX)/2
    averageY <- (minY + maxY)/2

    label_location <- rbind(label_location,c(Clusternumber, averageX, averageY))
  }
  label_location <- as.data.frame(label_location)
  colnames(label_location) <- c("Cluster", "Xpos", "Ypos")
  
  # use colourblind-friendly colours
  cluster_colours <- dittoColors()[seq_len(number_of_clusters)]

  q <- ggplot(cells_in_clusters, aes(x=Cell.X.Position, y=Cell.Y.Position))
  q <- q + geom_point(aes(color = Cluster))#, size = 0.01)
  q <- q + geom_text(data = label_location, aes(x = Xpos, y = Ypos, label = Cluster))
  q <- q + scale_color_manual(values=cluster_colours)
  q <- q + geom_point(data = cells_not_in_clusters,  colour = "black")#,size = 0.01)
  q <- q + xlab("Cell.X.Position") + ylab("Cell.Y.Position") +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")
  print(q)

  formatted_data_with_clusters <- formatted_data
  formatted_data_with_clusters$Cluster <- paste0("Cluster_", formatted_data_with_clusters$Cluster)
  formatted_data_with_clusters$Cluster[formatted_data_with_clusters$Cluster == "Cluster_NA"] <- "Free_cell"

  return(formatted_data_with_clusters)
}
