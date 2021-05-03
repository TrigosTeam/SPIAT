#' image_splitter
#'
#' @description Takes in a singlecellexperiment object from format_image_to_sce,
#' splits the image into specified sections and plots the different combinations
#' of markers within the image segments.
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param number_of_splits Number specifying the number of segments (e.g. 2 = 2x2, 3 = 3x3)
#' @param plot Specifies whether the splitted images should be printed in a pdf
#' @param cut_labels - specifies whether to plot where the image had been segmented
#' @param colour_vector List of colors. If specified, the colours will be used for plotting.
#' If NULL, colors will be generated automatically
#' @param x_position_min Integer used to specify the minimum x boundary to be splitted
#' @param x_position_max Integer used to specify the maximum x boundary to be splitted
#' @param y_position_min Integer used to specify the minimum y boundary to be splitted
#' @param y_position_max Integer used to specify the maximum y boundary to be splitted
#' @param column String specifying which column the colouring should be based on
#' @importFrom RColorBrewer brewer.pal
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#' @importFrom SummarizedExperiment colData
#' @importFrom stats complete.cases setNames
#' @importFrom dittoSeq dittoColors
#' @return A list of data.frames is returned
#' @examples
#' split_image <- image_splitter(SPIAT::formatted_image, number_of_splits=3, plot = FALSE)
#' @export

image_splitter <- function(sce_object, number_of_splits, plot = FALSE, cut_labels = TRUE, colour_vector = NULL,
                           x_position_min = NULL, x_position_max = NULL, y_position_min = NULL, y_position_max = NULL, 
                           column = "Cell.Type"){
    
    # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
    Cell.X.Position <- Cell.Y.Position <- NULL

    #turn the sce object name into string as a filename
    image_filename <- deparse(substitute(sce_object))

    #Reads the image file and deletes cell rows with NA positions
    cell_loc <- data.frame(colData(sce_object))
    cell_loc <- cell_loc[complete.cases(cell_loc),]
    
    #CHECK
    if (nrow(cell_loc) == 0) {
        stop("There are no cells in the data")
    }

    #Selects x and y region to plot
    if(is.null(x_position_min)){
        minX <- min(cell_loc$Cell.X.Position, na.rm = TRUE)
    }else{
        minX <- x_position_min
    }
    if(is.null(x_position_max)){
        maxX <- max(cell_loc$Cell.X.Position, na.rm = TRUE)
    }else{
        maxX <- x_position_max
    }
    if(is.null(y_position_min)){
        minY <- min(cell_loc$Cell.Y.Position, na.rm = TRUE)
    }else{
        minY <- y_position_min
    }
    if(is.null(y_position_max)){
        maxY <- max(cell_loc$Cell.Y.Position, na.rm = TRUE)
    }else{
        maxY <- y_position_max
    }

    #Partitions overall image into the manually selected x and y region
    manual_full_image <- cell_loc[minX < cell_loc$Cell.X.Position & cell_loc$Cell.X.Position <= maxX
                                  & minY < cell_loc$Cell.Y.Position & cell_loc$Cell.Y.Position <= maxY, ]
    if(isTRUE(plot)){
        number_markers <- length(unique(cell_loc[,column]))
        #Assigns colours to cell types based on user preference
        if(!is.null(colour_vector)){
            point_colours <- colour_vector
        }else{
        #Assigns colourblind-friendly colours    
            point_colours <- dittoColors()[seq_len(number_markers)]
        }

        #Plots partitioned full image
        full_image <- ggplot(manual_full_image, aes(x = Cell.X.Position, y = Cell.Y.Position)) +
            geom_point(aes(color = manual_full_image[,column]), size = 0.95) +
            scale_color_manual(values = point_colours) +
            guides(colour = guide_legend(title = "Cell Type", override.aes = list(size=1.0))) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "white"),
                  axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())

        #Plots border lines
        full_image <- full_image + geom_hline(aes_(yintercept = minY)) + geom_hline(aes_(yintercept = maxY)) +
            geom_vline(aes_(xintercept = minX)) + geom_vline(aes_(xintercept = maxX))

        #Labels border x and y positions of the full image
        if (isTRUE(cut_labels)){
            full_image <- full_image + geom_text(aes_string(x = 100, y = minY, label = minY), size = 3) +
                geom_text(aes_string(x = 100, y = maxY, label = maxY), size = 3) +
                geom_text(aes_string(x = minX, y = 50, label = minX), size = 3) +
                geom_text(aes_string(x = maxX, y = 50, label = maxX), size = 3)
        }

        pdf(paste(image_filename,"_image_split",".pdf",sep =""))
    }

    #Splits the range of x and y coordinates into n + 1 evenly spaced out lengths
    x_split <- seq(minX, maxX, length.out = number_of_splits + 1)
    y_split <- seq(minY, maxY, length.out = number_of_splits + 1)

    divided_image_obj <- list()
    for(y in seq_len(number_of_splits)){

        local_coor_y <- y_split[c(y+1, y)]

        #Rounds coordinates of cuts to nearest whole numberfor labelling purposes
        rounded_coor_y <- round(local_coor_y)

        #Creates horizontal gridlines on the original image file of where cuts are occurring
        if(isTRUE(plot) & y < number_of_splits) {

            full_image <- full_image + geom_hline(aes_(yintercept = max(local_coor_y)))

            #Labels cuts made along the y-axis
            if (isTRUE(cut_labels)) {
                full_image <- full_image + geom_text(aes_string(x = 100, y = max(local_coor_y), label = max(rounded_coor_y)), size = 3)
            }
        }

        temp_cell_loc <- cell_loc

        for(x in seq_len(number_of_splits)){

            local_coor_x <- x_split[c(x+1, x)]

            #Rounds coordinates of cuts to nearest whole numberfor labelling purposes
            rounded_coor_x <- round(local_coor_x)

            temp_cell_loc <- cell_loc

            #Creates vertical gridlines on the original image file of where cuts are occurring
            if(isTRUE(plot) & x > 1){

                full_image <- full_image + geom_vline(aes_(xintercept = min(local_coor_x)))

                #Labels cuts made along the x-axis
                if (isTRUE(cut_labels)){
                    full_image <- full_image + geom_text(aes_string(x = min(local_coor_x), y = 50, label = min(rounded_coor_x)), size = 3)
                }
            }
            #Locates all points in the data file fitting the given min/max parameters to be plotted i.e. splitting up the image.
            divided_image <- temp_cell_loc[min(local_coor_x) < temp_cell_loc$Cell.X.Position & temp_cell_loc$Cell.X.Position <= max(local_coor_x)
                                           & min(local_coor_y) < temp_cell_loc$Cell.Y.Position & temp_cell_loc$Cell.Y.Position <= max(local_coor_y), ]
            if(isTRUE(plot)){
                split_plot <- ggplot(divided_image, aes(x = Cell.X.Position, y = Cell.Y.Position)) +
                    geom_point(aes(colour = divided_image[,column]), size = 0.95) +
                    scale_colour_manual(values = point_colours) +
                    guides(colour = guide_legend(title = "Cell Type", override.aes = list(size=1.0))) +
                    ggtitle(paste("(", x, ", ", y, ")", sep = "")) +
                    theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_rect(fill = "white"),
                          axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank())

                print(split_plot)
            }

            divided_image_obj[[paste(image_filename,"r", x,"c", y, sep = "")]] <- divided_image
        }
    }
    if(isTRUE(plot)){
        print(full_image)
        dev.off()
        print("PDF saved successfully")
    }
    return(divided_image_obj)
}
