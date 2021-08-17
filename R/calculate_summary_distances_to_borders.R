#' calculate_summary_distances_to_borders
#'
#' @description Returns the mean, median and stardard deviation of the distances between a specified cell type to the borders
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param types_of_interest Vector of cell types to consider
#' @param column String specifying which column the interested cell types are from 
#' @importFrom stats median sd
#' @return A data.frame is returned
#' @export

calculate_summary_distances_to_borders <- function(sce_object, types_of_interest = NULL, column = NULL) {
  
  # CHECK if "Region" and "Distance.To.Border" columns exist
  if (is.null(sce_object$Region)){
    stop("Find the bordering cells first! Use the function identifying_bordering_cells_interactive()")
  }
  if (is.null(sce_object$Distance.To.Border)){
    stop("Find the min distances to the bordering cells first! Use the function calculate_min_distances_to_borders()")
  }
  
  # CHECK input
  if (is.null(types_of_interest)){
    stop("Please indicate the celltypes!")
  }
  if (is.null(column)){
    stop("Please indicate the columns of the interest!")
  }
  
  data <- data.frame(colData(sce_object))
  ##### data in #####
  data_of_interest_in <- dplyr::intersect(data[data[,column] %in% types_of_interest,],
                                          data[which(data$Region == "Inside"),])
  
  df.cols <- c("Min", "Max", "Mean",
               "Median", "St.dev")
  df <- data.frame(matrix(ncol=5,nrow=1, dimnames=list(NULL, df.cols)), stringsAsFactors = FALSE)
  
  min_d <- min(data_of_interest_in$Distance.To.Border)
  max_d <- max(data_of_interest_in$Distance.To.Border)
  mean_d <- mean(data_of_interest_in$Distance.To.Border)
  median_d <- median(data_of_interest_in$Distance.To.Border)
  st.dev_d <- sd(data_of_interest_in$Distance.To.Border)
  
  df <-  rbind(df[ ,df.cols], c(min_d, max_d, mean_d, median_d, st.dev_d))
  ##### data out #####
  data_of_interest_out <- dplyr::intersect(data[data[,column] %in% types_of_interest,],
                                           data[which(data$Region == "Outside"),])
  
  min_d <- min(data_of_interest_out$Distance.To.Border)
  max_d <- max(data_of_interest_out$Distance.To.Border)
  mean_d <- mean(data_of_interest_out$Distance.To.Border)
  median_d <- median(data_of_interest_out$Distance.To.Border)
  st.dev_d <- sd(data_of_interest_out$Distance.To.Border)
  
  
  df <-  rbind(df[ ,df.cols], c(min_d, max_d, mean_d, median_d, st.dev_d))
  rownames(df) <- c("empty","Inside","Outside")
  df <- df[-1,]   
  
  
  return(df)
}