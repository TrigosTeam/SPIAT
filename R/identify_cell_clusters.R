#' identify_cell_clusters
#'
#' @description Uses Euclidean distances to identify clusters of cells within a specified radius.
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param cell_types_of_interest Vector of phenotypes to consider
#' @param radius Integer specifying the radius of search.
#' @param feature_colname Column from which the cell types are selected
#' @param no_pheno Cell type corresponding to cells without a known phenotype (e.g. "None", "Other")
#' @import dplyr
#' @import ggplot2
#' @return A data.frame and a plot is returned
#' @examples
#' @export

# imported ggplo2 as interdependency of functions

identify_cell_clusters <- function(sce_object, method = "hierarchical", 
                                   cell_types_of_interest, radius, 
                                   min_cluster_size = 10,
                                   k=100,
                                   feature_colname, no_pheno = NULL) {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Cell.X.Position <- Cell.Y.Position <- Cluster <- Xpos <- Ypos <- NULL
  formatted_data <- get_colData(sce_object)
  
  ######remove cells without a phenotype
  if(!is.null(no_pheno)){
    formatted_data <- formatted_data[formatted_data[,feature_colname] != no_pheno, ]
  }
  
  #select cells to include if phenotypes of interest are specified
  if (!is.null(cell_types_of_interest)) {
    formatted_data <- formatted_data[formatted_data[,feature_colname] %in% cell_types_of_interest,]
  }
  
  #CHECK
  if (nrow(formatted_data) == 0) {
    stop("There are no cells in data/no cells for the phenotypes of interest")
  }
  
  
  # hierarchical clustering
  if (method == "hierarchical"){
    rownames(formatted_data) <- formatted_data$Cell.ID
    sim_close <- - apcluster::negDistMat(formatted_data[,c("Cell.X.Position", "Cell.Y.Position")])
    
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
        h <- stats::hclust(stats::as.dist(sim_close), method="single")
        
        local_clusters <- stats::cutree(h, h = 0.5)
        
        formatted_data$Cluster <- as.character(local_clusters[match(formatted_data$Cell.ID, names(local_clusters))])
        
      } else {
        formatted_data$Cluster <- NA
        #stop("The radius specified may be too small, no clusters were found")
      }
    } else {
      formatted_data$Cluster <- NA  
      #stop("The radius specified may be too small, no clusters were found")
    }
    
  }
  else if (method == "dbscan"){
    cell_cords <- formatted_data[,c("Cell.X.Position", "Cell.Y.Position")]
    #Use dbscan to generate clusters
    db <- dbscan::dbscan(cell_cords, eps = radius, minPts = min_cluster_size)
    #since dbscan outputs cluster 0 as noise, we add 1 to all cluster numbers to keep it consistent
    formatted_data$Cluster <- factor(db$cluster + 1)
  }
  else if (method == "rphenograph"){
    cell_cords <- formatted_data[,c("Cell.X.Position", "Cell.Y.Position")]
    Rphenograph_out <- Rphenograph(cell_cords, k = k)
    formatted_data$Cluster <- factor(membership(Rphenograph_out[[2]]))
  }
  else {
    stop("Please select a valid clustering method, current options: dbscan")
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
  cluster_colours <- dittoSeq::dittoColors()[seq_len(number_of_clusters)]
  
  q <- ggplot(cells_in_clusters, aes(x=Cell.X.Position, y=Cell.Y.Position))
  q <- q + geom_point(aes(color = Cluster))#, size = 0.01)
  q <- q + geom_text(data = label_location, aes(x = Xpos, y = Ypos, label = Cluster))
  q <- q + scale_color_manual(values=cluster_colours)
  if(dim(cells_not_in_clusters)[1]!=0){
    q <- q + geom_point(data = cells_not_in_clusters,  colour = "black")#,size = 0.01)
  }
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
