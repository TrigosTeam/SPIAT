#' calculate_summary_distances_of_cells_to_borders
#'
#' @description Returns the mean, median and standard deviation of the distances
#'   between a specified cell type to the border.
#' @param spe_object SpatialExperiment object. It should contain information of
#'   tumour structure and cell distances to tumour border (`colData()` has
#'   `Region` and `Distance.To.Border` columns).
#' @param cell_types_of_interest String Vector of cell types to consider.
#' @param feature_colname String specifying which column the interested cell
#'   types are from.
#' @return A data.frame is returned
#' @export
#' @examples
#' spe_border <- identify_bordering_cells(SPIAT::defined_image,
#' reference_cell = "Tumour", feature_colname = "Cell.Type", n_to_exclude = 10)
#' spe_dist <- calculate_distance_to_tumour_margin(spe_border)
#' spe_structure <- define_structure(spe_dist, names_of_immune_cells =
#' c("Immune1","Immune2","Immune3"), feature_colname = "Cell.Type",
#' n_margin_layers = 5)
#' calculate_summary_distances_of_cells_to_borders(spe_structure,
#' cell_types_of_interest = c("Immune1","Immune3"),feature_colname = "Cell.Type")

calculate_summary_distances_of_cells_to_borders <- function(
        spe_object,  cell_types_of_interest,  feature_colname = "Cell.Type") {
    
    # CHECK if "Region" and "Distance.To.Border" columns exist
    if (is.null(spe_object$Region)){
        stop("Find the bordering cells first! Use the function identifying_bordering_cells_interactive()")
    }
    if (is.null(spe_object$Distance.To.Border)){
        stop("Find the min distances to the bordering cells first! Use the function calculate_min_distances_to_borders()")
    }
    
    # CHECK input
    if (is.null(cell_types_of_interest)){
        stop("Please indicate the cell types!")
    }
    if (is.null(feature_colname)){
        stop("Please indicate the column name of the cell types!")
    }
    
    data <- data.frame(SummarizedExperiment::colData(spe_object))
    
    # define a function to get the statistics of the distances
    summarise_dist <- function(data){
        if (dim(data)[1] == 0){
            min_d <- max_d <- mean_d <- median_d <- st.dev_d <- NA
        } 
        else {
            min_d <- min(data$Distance.To.Border, na.rm = TRUE)
            max_d <- max(data$Distance.To.Border, na.rm = TRUE)
            mean_d <- mean(data$Distance.To.Border, na.rm = TRUE)
            median_d <- stats::median(data$Distance.To.Border, na.rm = TRUE)
            st.dev_d <- stats::sd(data$Distance.To.Border, na.rm = TRUE)
        }
        return(c(min_d, max_d, mean_d, median_d, st.dev_d))
    }
    ##### data in #####
    data_of_interest_in <- 
        data[which((data[[feature_colname]] %in% cell_types_of_interest) 
                   & (data$Region == "Inside")),]
    
    df.cols <- c("Cell.Type", "Location", "Min", "Max", "Mean",
                 "Median", "St.dev")
    df <- vector()
    
    sum_d <- summarise_dist(data_of_interest_in)
    
    df <-  rbind(df, c("All_cell_types_of_interest", "Tumor_area", sum_d))
    ##### data out #####
    data_of_interest_out <- 
        data[which((data[[feature_colname]] %in% cell_types_of_interest) 
                   & (data$Region == "Outside")),]
    
    sum_d <- summarise_dist(data_of_interest_out)
    
    df <-  rbind(df, c("All_cell_types_of_interest", "Stroma", sum_d))
    
    for(type in cell_types_of_interest){
        data_of_interest_in <- data[which((data[[feature_colname]] %in% type) &
                                              (data$Region == "Inside")),]
        sum_d <- summarise_dist(data_of_interest_in)
        
        df <-  rbind(df, c(type, "Tumor_area", sum_d))
        ##### data out #####
        data_of_interest_out <- data[which((data[[feature_colname]] %in% type) &
                                               (data$Region == "Outside")),]
        sum_d <- summarise_dist(data_of_interest_out)
        
        df <-  rbind(df, c(type, "Stroma", sum_d))
    }
    colnames(df) <- c("Cell.Type", "Area", "Min_d", "Max_d", "Mean_d", 
                      "Median_d", "St.dev_d")
    df <- as.data.frame(df)
    df$Min_d <- as.numeric(df$Min_d)
    df$Max_d <- as.numeric(df$Max_d)
    df$Mean_d <- as.numeric(df$Mean_d)
    df$Median_d <- as.numeric(df$Median_d)
    df$St.dev_d <- as.numeric(df$St.dev_d)
    return(df)
}
