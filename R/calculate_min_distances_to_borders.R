#' calculate_min_distance_to_borders
#'
#' @description Returns the sce_object with the minimum distance from immune cells to the identified bordering cells
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param split integer specifying the number of subimages to split; if NULL, do not split
#' @import dplyr
#' @import magrittr
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble rownames_to_column
#' @importFrom stats complete.cases
#' @importFrom spatstat.geom crossdist.default
#' @return An sce_object is returned
#' @export

calculate_min_distance_to_borders <- function(sce_object, split = NULL){
  
  #CHECK if the user has found the bordering cells yet
  if (is.null(sce_object$Region)){
    stop("Please find the bordering cells first!")
  }
  
  #Reads the image file and deletes cell rows with NA positions
  dat <- data.frame(colData(sce_object))
  dat <- dat[complete.cases(dat),]
  dat<- dat %>% rownames_to_column("Cell.ID") #convert rowname to column
  
  #CHECK
  if (nrow(dat) == 0) {
    stop("There are no cells")
  }
  
  dat <- dat[, c("Cell.ID","Region", "Cell.X.Position", "Cell.Y.Position")]
  dat <- dat[dat$Region != "",]
  
  first_dataset <- dat[which(dat$Region == "Border"),]
  second_dataset <- dat[which(dat$Region != "Border"),]
  
  if (is.null(split)){
    dist_matrix <- crossdist.default(first_dataset[,"Cell.X.Position"], first_dataset[,"Cell.Y.Position"], 
                                     second_dataset[,"Cell.X.Position"], second_dataset[,"Cell.Y.Position"] )
    colnames(dist_matrix)<-second_dataset$Cell.ID
    rownames(dist_matrix)<-first_dataset$Cell.ID
    
    
    border_cell_ids <- dat[which(dat$Region == "Border"),"Cell.ID"]
    border_cells <- rep(0, length(border_cell_ids))
    names(border_cells) <- border_cell_ids
    
    other_cells<- apply(dist_matrix,2,min)
    
    min_distance <- c(other_cells,border_cells)
    rownames(dat) <- dat$Cell.ID
    dat[names(min_distance),"Distance.To.Border"] <- min_distance
    colData(sce_object)$Distance.To.Border <- dat[,"Distance.To.Border"]
  }
  else{
    second_sce <- format_colData_to_sce(second_dataset)
    
    # split the data
    splits <- image_splitter(second_sce, split)
    border_cell_ids <- dat[which(dat$Region == "Border"),"Cell.ID"]
    border_cells <- rep(0, length(border_cell_ids))
    names(border_cells) <- border_cell_ids
    min_distance <- c(border_cells)
    for (sub_image in splits){
      dist_matrix <- crossdist.default(first_dataset[,"Cell.X.Position"], first_dataset[,"Cell.Y.Position"], 
                                       sub_image[,"Cell.X.Position"], sub_image[,"Cell.Y.Position"] )
      colnames(dist_matrix)<-sub_image$Cell.ID
      rownames(dist_matrix)<-first_dataset$Cell.ID
      
      other_cells<- apply(dist_matrix,2,min)
      min_distance <- c(min_distance, other_cells)
    }
    rownames(dat) <- dat$Cell.ID
    dat[names(min_distance),"Distance.To.Border"] <- min_distance
    colData(sce_object)$Distance.To.Border <- dat[,"Distance.To.Border"]
  }
  
  return(sce_object)
}