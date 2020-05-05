#' marker_surface_plot
#'
#' @description Generates a 3D surface plot of the level of the selected
#' marker. Note that the image is blurred based on the 'num_splits' parameter.
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param num_splits Integer specifying the number of splits on the image, higher
#' splits equal to higher resolution. Recommendation: 10-100
#' @param marker Marker to plot
#' @param x_position_min Integer specifying the minimum x boundary to be splitted
#' @param x_position_max Integer specifying the maximum x boundary to be splitted
#' @param y_position_min Integer specifying the minimum y boundary to be splitted
#' @param y_position_max Integer specifying the maximum y boundary to be splitted
#' @import dplyr
#' @import SingleCellExperiment
#' @importFrom tibble rownames_to_column
#' @importFrom stats aggregate
#' @importFrom plotly plot_ly add_surface
#' @export

# %>% operator is in package 'magrittr' but imported by dplyr
# colData() is in package 'SummarizedExperiment' but imported by SingleCellExperiment

marker_surface_plot <- function(sce_object, num_splits, marker, x_position_min = NULL, x_position_max = NULL,
                                y_position_min = NULL, y_position_max = NULL){

    formatted_data <- data.frame(colData(sce_object))

    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

    expression_matrix <- assay(sce_object)

    markers <- rownames(expression_matrix)
    
    #CHECK
    if (is.element(marker, markers) == FALSE) {
        stop("Data does not contain marker specified")
    }
    
    cell_ids <- colnames(expression_matrix)

    rownames(expression_matrix) <- NULL
    colnames(expression_matrix) <- NULL
    expression_matrix_t <- t(expression_matrix)
    expression_df <- data.frame(expression_matrix_t)
    colnames(expression_df) <- markers

    formatted_data <- cbind(formatted_data, expression_df)
    formatted_data <- formatted_data[complete.cases(formatted_data),]

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
    for (y in 1:num_splits){
        local_coor_y <- y_split[c(y+1, y)]
        #print(local_coor_y)

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

    for (x in 1:num_splits){
        local_coor_x <- x_split[c(x+1, x)]
        # print(local_coor_x)

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

    #create a df with only the expression level of a single marker of interest and the coordinates
    df <- aggregate(formatted_data[,marker], by=list(xcord=formatted_data$split.X, ycord=formatted_data$split.Y), FUN=mean)

    #initialize a matrix for surface plot, dim=num_splits^2
    my_matrix <- matrix(nrow = num_splits, ncol=num_splits)

    #populate matrix with values from df
    for (x in 1:num_splits){

        for (y in 1:num_splits){

            #select the row with the xcord and ycord
            row <- df[df[, "xcord"] == x & df[, "ycord"] == y, ]

            #if there is expression in that coordinate, assign it to matrix
            if (nrow(row) == 1) {
                my_matrix[x,y] <- row$x
            }
            else {
                my_matrix[x,y] <- 0
            }

        }
    }

    #somehow transposing matrix matches the plots from before, without transposing it is mirrored...
    my_matrix <- t(my_matrix)

    #rename the matrix
    mean_marker_level <- my_matrix

    p <- plot_ly(z = ~mean_marker_level) %>% add_surface()
    print(p)

}
