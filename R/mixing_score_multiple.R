#' mixing_score_multiple
#'
#' @description Produces a list of dataframes that contains all the mixing scores respective to each reference cell. Each dataframe is corresponding to one reference cell type 
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param reference_marker Markers to be used as the reference cells. 
#' @param target_marker Markers to be used as the target markers. The function will ignore instances where reference marker is the same as target marker. 
#' @param column String specifying the column where the selected cell types are from
#' @param radius (OPTIONAL) The maximum radius around a reference marker for another cell to be considered an interaction.
#' @import dplyr
#' @importFrom SummarizedExperiment colData 
#' @importFrom tibble rownames_to_column
#' @importFrom dbscan frNN
#' @return A list of dataframes
#' @examples
#' @export

mixing_score_multiple <- function(sce_object, reference_marker, target_marker, radius=20, column = "Phenotype")
{
  formatted_data <- data.frame(colData(sce_object))
  formatted_data <- formatted_data[complete.cases(formatted_data), 
  ]
  formatted_data <- formatted_data %>% rownames_to_column("Cell.ID")
  df.cols <- c("Reference", "Number_of_reference_cells",
               "Number_of_target_cells", "Reference_target_interaction",
               "Reference_reference_interaction", "Mixing_score", 
               "Normalised_mixing_score")
  df <- data.frame(matrix(ncol=7,nrow=1, dimnames=list(NULL, df.cols)), stringsAsFactors = FALSE)
  all.df <- list()
  for (i in reference_marker) {
    reference_cells <- formatted_data[which(formatted_data[,column] == i), ]
    reference_cell_cords <- reference_cells[, c("Cell.ID", "Cell.X.Position", 
                                                "Cell.Y.Position")]
    dataframe <- remove_rownames(reference_cell_cords)
    dataframe <- dataframe %>% column_to_rownames("Cell.ID")
    reference_cell_cords <- reference_cells[, c( "Cell.X.Position", 
                                                "Cell.Y.Position")]
    reference_reference_result <- frNN(reference_cell_cords, eps = radius, sort = FALSE)
    reference_reference_interactions <- sum(rapply(reference_reference_result$id, length))
    dataframe[,i] <- rapply(reference_reference_result$id, length)
    dataframe[which(dataframe[, i] == 0),i] <- NA
    dataframe[, "No.Target"] <- 0
    
    n_target <- 0
    for (j in target_marker) {
      
      if (i == j) {next}
      tryCatch({
        target_cells <- formatted_data[which( formatted_data[,column] == j), ]
        n_target <- n_target +nrow(target_cells)
        if (nrow(reference_cells) == 0) {
          print(paste("There are no unique reference cells of specified marker", i, "for target cell", j))
        }
        if (nrow(target_cells) == 0) {
          print(paste("There are no unique target cells of specified marker", j, "for reference cell", i))
        }
       
        target_cell_cords <- target_cells[, c("Cell.ID", "Cell.X.Position", 
                                              "Cell.Y.Position")]
      	target_cell_cords <- remove_rownames(target_cell_cords)
      	target_cell_cords <- target_cell_cords %>% column_to_rownames("Cell.ID")
      	reference_target_result <- frNN(target_cell_cords, eps = radius, 
      	                                query = reference_cell_cords, sort = FALSE)
      	dataframe[,j] <- rapply(reference_target_result$id, length)
      	dataframe[, "No.Target"] <- dataframe[,j]+dataframe[, "No.Target"] 

       }, error=function(e){})
    }

    dataframe[,"Mixing.Score"] <- dataframe[,"No.Target"]/dataframe[,i]
    mixing_score <- sum(dataframe[,"No.Target"])/reference_reference_interactions
    normalised_mixing_score <- mixing_score * nrow(dataframe) / nrow(target_cells)
    p <- density(as.numeric(dataframe[,"Mixing.Score"]),na.rm=T)
    df <-  rbind(df[ ,df.cols], 
                 c(i, nrow(dataframe), n_target, sum(dataframe[,"No.Target"]), 
                   reference_reference_interactions, mixing_score, normalised_mixing_score))
    all.df[[i]]$data <- dataframe
    all.df[[i]]$density <- p
  }
  df[,3:7] <- sapply(df[,3:7],as.numeric)
  
  
  df <- df[-1,]
  all.df[["Summary"]] <- df
  return(all.df)
}
