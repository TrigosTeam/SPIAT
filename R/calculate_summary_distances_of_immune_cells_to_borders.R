#' calculate_summary_distances_of_cells_to_borders
#'
#' @description Returns the mean, median and stardard deviation of the distances
#'   between a specified cell type to the borders.
#' @param sce_object SingleCellExperiment object in the form of the output of
#'   format_image_to_sce.
#' @param cell_types_of_interest String Vector of cell types to consider.
#' @param feature_colname String specifying which column the interested cell
#'   types are from.
#' @return A data.frame is returned
#' @export
#' @examples 
#' sce_border <- identify_bordering_cells(SPIAT::defined_image, reference_cell = "Tumour",
#' feature_colname = "Cell.Type", n_to_exclude = 10)
#' sce_dist <- calculate_distance_to_tumour_margin(sce_border)
#' sce_structure <- define_structure(sce_dist, names_of_immune_cells = c("Immune1","Immune2","Immune3"),
#' feature_colname = "Cell.Type", n_margin_layers = 5)
#' calculate_summary_distances_of_cells_to_borders(sce_structure, 
#' cell_types_of_interest = c("Immune1","Immune3"),feature_colname = "Cell.Type")

calculate_summary_distances_of_cells_to_borders <- function(sce_object, 
                                                            cell_types_of_interest, 
                                                            feature_colname = "Cell.Type") {
  
  # CHECK if "Region" and "Distance.To.Border" columns exist
  if (is.null(sce_object$Region)){
    stop("Find the bordering cells first! Use the function identifying_bordering_cells_interactive()")
  }
  if (is.null(sce_object$Distance.To.Border)){
    stop("Find the min distances to the bordering cells first! Use the function calculate_min_distances_to_borders()")
  }
  
  # CHECK input
  if (is.null(cell_types_of_interest)){
    stop("Please indicate the celltypes!")
  }
  if (is.null(feature_colname)){
    stop("Please indicate the columns of the interest!")
  }
  
  data <- data.frame(colData(sce_object))
  ##### data in #####
  data_of_interest_in <- data[which((data[[feature_colname]] %in% cell_types_of_interest) 
                                    & (data$Region == "Inside")),]
  
  df.cols <- c("Cell.Type", "Location", "Min", "Max", "Mean",
               "Median", "St.dev")
  df <- vector()
  
  min_d <- min(data_of_interest_in$Distance.To.Border)
  max_d <- max(data_of_interest_in$Distance.To.Border)
  mean_d <- mean(data_of_interest_in$Distance.To.Border)
  median_d <- stats::median(data_of_interest_in$Distance.To.Border)
  st.dev_d <- stats::sd(data_of_interest_in$Distance.To.Border)
  
  df <-  rbind(df, c("All_cell_types_of_interest", "Tumor_area", min_d, max_d, mean_d, median_d, st.dev_d))
  ##### data out #####
  data_of_interest_out <- data[which((data[[feature_colname]] %in% cell_types_of_interest) 
                                     & (data$Region == "Outside")),]
  
  min_d <- min(data_of_interest_out$Distance.To.Border)
  max_d <- max(data_of_interest_out$Distance.To.Border)
  mean_d <- mean(data_of_interest_out$Distance.To.Border)
  median_d <- stats::median(data_of_interest_out$Distance.To.Border)
  st.dev_d <- stats::sd(data_of_interest_out$Distance.To.Border)
  
  df <-  rbind(df, c("All_cell_types_of_interest", "Stroma", min_d, max_d, mean_d, median_d, st.dev_d))
  
  for(type in cell_types_of_interest){
    data_of_interest_in <- data[which((data[[feature_colname]] %in% type) & (data$Region == "Inside")),]
    min_d <- min(data_of_interest_in$Distance.To.Border)
    max_d <- max(data_of_interest_in$Distance.To.Border)
    mean_d <- mean(data_of_interest_in$Distance.To.Border)
    median_d <- median(data_of_interest_in$Distance.To.Border)
    st.dev_d <- sd(data_of_interest_in$Distance.To.Border)
    
    df <-  rbind(df, c(type, "Tumor_area", min_d, max_d, mean_d, median_d, st.dev_d))
    ##### data out #####
    data_of_interest_out <- data[which((data[[feature_colname]] %in% type) & (data$Region == "Outside")),]
    
    min_d <- min(data_of_interest_out$Distance.To.Border)
    max_d <- max(data_of_interest_out$Distance.To.Border)
    mean_d <- mean(data_of_interest_out$Distance.To.Border)
    median_d <- stats::median(data_of_interest_out$Distance.To.Border)
    st.dev_d <- stats::sd(data_of_interest_out$Distance.To.Border)
    
    df <-  rbind(df, c(type, "Stroma", min_d, max_d, mean_d, median_d, st.dev_d))
  }
  colnames(df) <- c("Cell.Type", "Area", "Min_d", "Max_d", "Mean_d", "Median_d", "St.dev_d")
  df <- as.data.frame(df)
  df$Min_d <- as.numeric(df$Min_d)
  df$Max_d <- as.numeric(df$Max_d)
  df$Mean_d <- as.numeric(df$Mean_d)
  df$Median_d <- as.numeric(df$Median_d)
  df$St.dev_d <- as.numeric(df$St.dev_d)
  return(df)
}
