#' plot_marker_level_heatmap
#'
#' @description Blurs the image by splitting the images into small squares. The
#'   marker levels are then averaged within each square. All cells are
#'   considered, regardless of phenotype status.
#'
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param marker String. Marker to plot.
#' @param num_splits Integer specifying the blurring level (number of splits)
#'   for the image. Higher numbers result in higher resolution.
#' @import dplyr
#' @import ggplot2
#' @return A plot is returned
#' @examples
#' plot_marker_level_heatmap(SPIAT::simulated_image, num_splits = 100, "Tumour_marker")
#' @export

plot_marker_level_heatmap <- function(spe_object, num_splits, marker){
    
    # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
    xcord <- ycord <- NULL

    formatted_data <- get_colData(spe_object)

    intensity_matrix <- SummarizedExperiment::assay(spe_object)

    markers <- rownames(intensity_matrix)
    
    #CHECK
    if (is.element(marker, markers) == FALSE) {
        stop("The marker specified is not in the data")
    }
    
    cell_ids <- colnames(intensity_matrix)

    rownames(intensity_matrix) <- NULL
    colnames(intensity_matrix) <- NULL
    intensity_matrix_t <- t(intensity_matrix)
    intensity_df <- data.frame(intensity_matrix_t)
    colnames(intensity_df) <- markers

    formatted_data <- cbind(formatted_data, intensity_df)
    formatted_data <- formatted_data[stats::complete.cases(formatted_data),]

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

    heatmap_title <- paste(marker, "level")

    #create a df with only the intensity level of a single marker of interest and the coordinates
    df <- stats::aggregate(formatted_data[,marker], by=list(xcord=formatted_data$split.X, ycord=formatted_data$split.Y), FUN=mean)

    p <- ggplot(df, aes(xcord, ycord, fill=x)) + geom_tile()
    p <- p + scale_fill_gradient(low="white", high="red")
    p <- p + xlab("x position") + ylab("y position")
    p <- p + labs(fill = "Mean intensity level") + ggtitle(heatmap_title)
    p <- p + theme(panel.background = element_rect(fill = "grey", colour = "grey", linetype = "solid"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
    methods::show(p)
    return(p)
}
