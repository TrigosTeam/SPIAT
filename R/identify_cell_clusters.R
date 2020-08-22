#' identify_cell_clusters
#'
#' @description Uses Euclidean distances to identify clusters of cells within a specified radius.
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param phenotypes_of_interest Vector of phenotypes to consider
#' @param radius Integer specifying the radius of search.
#' @import dplyr
#' @importFrom SummarizedExperiment colData assay
#' @importFrom tibble rownames_to_column
#' @importFrom stats complete.cases hclust cutree as.dist
#' @importFrom apcluster negDistMat
#' @importFrom dittoSeq dittoColors
#' @import ggplot2
#' @export

# imported ggplo2 as interdependency of functions

identify_cell_clusters <- function(sce_object, phenotypes_of_interest, radius) {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Cell.X.Position <- Cell.Y.Position <- Cluster <- Xpos <- Ypos <- NULL
  
  formatted_data <- data.frame(colData(sce_object))
  formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

  expression_matrix <- assay(sce_object)

  markers <- rownames(expression_matrix)
  cell_ids <- colnames(expression_matrix)

  rownames(expression_matrix) <- NULL
  colnames(expression_matrix) <- NULL
  expression_matrix_t <- t(expression_matrix)
  expression_df <- data.frame(expression_matrix_t)
  colnames(expression_df) <- markers

  formatted_data <- cbind(formatted_data, expression_df)
  formatted_data <- formatted_data[complete.cases(formatted_data),]

  ######remove cells without a phenotype
  formatted_data <- formatted_data[formatted_data$Phenotype != "OTHER", ]
  #cell_cords <- formatted_data[,c("Cell.X.Position", "Cell.Y.Position")]

  #select cells to include if phenotypes of interest are specified
  if (!is.null(phenotypes_of_interest)) {
    formatted_data <- formatted_data[formatted_data$Phenotype %in% phenotypes_of_interest,]
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
        stop("The radius specified may be too small, no clusters were found")
      }
  } else {
      stop("The radius specified may be too small, no clusters were found")
    }

  #get cells assigned to clusters
  cells_in_clusters <- formatted_data[complete.cases(formatted_data),]
  
  #get number_of_clusters
  number_of_clusters <- length(unique(cells_in_clusters$Cluster))

  #label the Cluster centre by averaging x and y

  label_location <- vector()
  for (Clusternumber in 1:number_of_clusters) {
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
  colours <- dittoColors()[1:number_of_clusters]

  q <- ggplot(cells_in_clusters, aes(x=Cell.X.Position, y=Cell.Y.Position))
  q <- q + geom_point(aes(color = Cluster), size = 0.01)
  q <- q + geom_text(data = label_location, aes(x = Xpos, y = Ypos, label = Cluster))
  q <- q + xlab("Cell.X.Position") + ylab("Cell.Y.Position")
  q <- q + scale_color_manual(values=colours) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")
  print(q)

  formatted_data_with_clusters <- formatted_data
  formatted_data_with_clusters$Cluster <- paste0("Cluster_", formatted_data_with_clusters$Cluster)

  return(formatted_data_with_clusters)
}
