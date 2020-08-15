#' plot_marker_level_heatmap
#'
#' @description Blurres the image by splitting the images into small squares.
#' The marker levels are then averaged within each square. All cells are considered,
#' regardless of phenotype status
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param marker Marker to plot
#' @param num_splits Integer specifying the blurring level (number of splits) for the image.
#' Higher numbers result in higher resolution.
#' @import dplyr
#' @import ggplot2
#' @importFrom SummarizedExperiment colData assay
#' @importFrom tibble rownames_to_column
#' @importFrom stats aggregate
#' @export

plot_marker_level_heatmap <- function(sce_object, num_splits, marker){

    formatted_data <- data.frame(colData(sce_object))

    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

    expression_matrix <- assay(sce_object)

    markers <- rownames(expression_matrix)
    
    #CHECK
    if (is.element(marker, markers) == FALSE) {
        stop("The marker specified is not in the data")
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

    minX <- min(formatted_data$Cell.X.Position, na.rm = TRUE)
    maxX <- max(formatted_data$Cell.X.Position, na.rm = TRUE)
    minY <- min(formatted_data$Cell.Y.Position, na.rm = TRUE)
    maxY <- max(formatted_data$Cell.Y.Position, na.rm = TRUE)

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

    heatmap_title <- paste(marker, "level")

    #create a df with only the expression level of a single marker of interest and the coordinates
    df <- aggregate(formatted_data[,marker], by=list(xcord=formatted_data$split.X, ycord=formatted_data$split.Y), FUN=mean)

    p <- ggplot(df, aes(xcord, ycord, fill=x)) + geom_tile()
    p <- p + scale_fill_gradient(low="white", high="red")
    p <- p + xlab("x position") + ylab("y position")
    p <- p + labs(fill = "Mean expression level") + ggtitle(heatmap_title)
    p <- p + theme(panel.background = element_rect(fill = "grey", colour = "grey", linetype = "solid"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
    print(p)
}
