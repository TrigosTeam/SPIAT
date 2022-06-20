#' marker_surface_plot_stack
#'
#' @description Generates stacked 3D surface plots showing normalized intensity
#'   level of specified markers.
#'
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param num_splits Integer specifying the number of splits on the image,
#'   higher splits equal to higher resolution. Recommendation: 10-100.
#' @param markers_to_plot Vector of marker names for plotting.
#' @param sep Integer specifying the distance separation between each surface
#'   plot. We recommend values in the 1-2 range.
#' @param x_position_min Integer specifying the minimum x boundary to be
#'   splitted.
#' @param x_position_max Integer specifying the maximum x boundary to be
#'   splitted.
#' @param y_position_min Integer specifying the minimum y boundary to be
#'   splitted.
#' @param y_position_max Integer specifying the maximum y boundary to be
#'   splitted.
#' @import dplyr
#' @return A plot is returned
#' @examples
#' marker_surface_plot_stack(SPIAT::simulated_image, num_splits=15, 
#' markers=c("Tumour_marker", "Immune_marker4"))
#' @export

marker_surface_plot_stack <- function(spe_object, num_splits, markers_to_plot, sep = 1,
                                      x_position_min = NULL, x_position_max = NULL,
                                      y_position_min = NULL, y_position_max = NULL){
    
    #CHECK
    intensity_df <- SummarizedExperiment::assay(spe_object)
    markers <- rownames(intensity_df)
    if (!all(markers_to_plot %in% markers)) {
        stop("One or more markers specified cannot be found")
    }
    
    # format data
    formatted_data <- bind_info(spe_object)
    formatted_data$split.X <- 0
    formatted_data$split.Y <- 0
    
    
    #Selects x and y region to plot
    if(is.null(x_position_min)){
        minX <- min(formatted_data$Cell.X.Position, na.rm = TRUE)
    }else{
        minX <- x_position_min
    }
    if(is.null(x_position_max)){
        maxX <- max(formatted_data$Cell.X.Position, na.rm = TRUE)
    }else{
        maxX <- x_position_max
    }
    if(is.null(y_position_min)){
        minY <- min(formatted_data$Cell.Y.Position, na.rm = TRUE)
    }else{
        minY <- y_position_min
    }
    if(is.null(y_position_max)){
        maxY <- max(formatted_data$Cell.Y.Position, na.rm = TRUE)
    }else{
        maxY <- y_position_max
    }
    
    #Splits the range of x and y coordinates
    #into n + 1 evenly spaced out lengths
    x_split <- seq(minX, maxX, length.out = num_splits + 1)
    y_split <- seq(minY, maxY, length.out = num_splits + 1)
    
    #Creates matrix of the locations of x and y cuts to the image
    split_occurrence <- cbind(x_split, y_split)
    
    #obtain the x and y coordinates on a heatmap for every cell based on number of splits
    for (y in seq_len(num_splits)){
        local_coor_y <- y_split[c(y+1, y)]
        
        #grab the cells in the range
        result <- formatted_data[min(local_coor_y) < formatted_data$Cell.Y.Position & formatted_data$Cell.Y.Position <= max(local_coor_y), ]
        if(y == 1){
            extra_row <- formatted_data[formatted_data$Cell.Y.Position == min(local_coor_y), ]
            result <- rbind(result, extra_row)
        }
        
        if(nrow(result) > 0) {
            result$split.Y <- y
            formatted_data[match(result$Cell.ID,formatted_data$Cell.ID),] <- result
        }
    }
    
    for (x in seq_len(num_splits)){
        local_coor_x <- x_split[c(x+1, x)]
        
        #grab the cells in the range
        result <- formatted_data[min(local_coor_x) < formatted_data$Cell.X.Position & formatted_data$Cell.X.Position <= max(local_coor_x), ]
        if(x == 1){
            extra_row <- formatted_data[formatted_data$Cell.X.Position == min(local_coor_x), ]
            result <- rbind(result, extra_row)
        }
        
        if(nrow(result) > 0) {
            result$split.X <- x
            formatted_data[match(result$Cell.ID,formatted_data$Cell.ID),] <- result
        }
    }
    
    #start plotting the surface plots
    p <- plotly::plot_ly()
    
    #value to separate the plots
    i <- 0
    
    for (marker in markers_to_plot) {
        #skip DAPI intensities
        if (marker == "DAPI"){
            next
        }
        
        #create a df with only the intensity level of a single marker of interest and the coordinates
        df <- stats::aggregate(formatted_data[,marker], by=list(xcord=formatted_data$split.X, ycord=formatted_data$split.Y), FUN=mean)
        
        #initialize a matrix for surface plot, dim=num_splits^2
        my_matrix <- matrix(nrow = num_splits, ncol=num_splits)
        
        #populate matrix with values from df
        for (x in seq_len(num_splits)){
            
            for (y in seq_len(num_splits)){
                
                #select the row with the xcord and ycord
                row <- df[df[, "xcord"] == x & df[, "ycord"] == y, ]
                
                #if there is intensity in that coordinate, assign it to matrix
                if (nrow(row) == 1) {
                    my_matrix[x,y] <- row$x
                }
                else {
                    my_matrix[x,y] <- 0
                }
            }
        }
        
        #function to scale values between 0.1-1.1
        normalize <- function(x){
            return((x-min(x, na.rm=TRUE))/(max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
        }
        
        #scale the matrix
        my_matrix <- normalize(my_matrix)
        
        #Transposing so resulting plot matches that of the other plotting functions
        my_matrix <- t(my_matrix)
        
        #add the value to separate plots
        my_matrix <- my_matrix + i
        
        # use colourblind-friendly colour
        hexcol <- dittoSeq::dittoColors()[i + 1]
        rgbcol_mat <- grDevices::col2rgb(hexcol)
        rgbcol <- paste0(rgbcol_mat[,1], collapse=",")
        rgbcol <- paste0("rgb(", rgbcol, ")") 
        
        #add the surface
        p <- plotly::add_trace(p, z = my_matrix, type = "surface", 
                               colorscale = list(c(0,1),c("rgb(255,255,255)", rgbcol)), 
                               colorbar=list(title=marker))
        i <- i + sep
    }
    p
}
