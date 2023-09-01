#' plot_cell_categories
#'
#' @description Produces a scatter plot of the cells of their x-y positions in
#'   the tissue. Cells are coloured categorically by phenotype. Cells not part
#'   of the phenotypes of interest will be coloured "lightgrey".
#'
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param categories_of_interest Vector of cell categories to be coloured.
#' @param colour_vector Vector specifying the colours of each cell phenotype.
#' @param feature_colname String specifying the column the cell categories
#'   belong to.
#' @param cex Numeric. The size of the plot points. Default is 1.
#' @param layered Boolean. Whether to plot the cells layer by layer (cell
#'   categories). By default is FALSE.
#' @import dplyr
#' @import ggplot2
#' @importFrom rlang .data
#' @return A plot is returned
#' @examples
#' categories_of_interest <- c("Tumour", "Immune1","Immune2","Immune3")
#' colour_vector <- c("red","darkblue","blue","darkgreen")
#' plot_cell_categories(SPIAT::defined_image, categories_of_interest, colour_vector,
#' feature_colname = "Cell.Type")
#' @export

plot_cell_categories <- function(spe_object, categories_of_interest = NULL, 
                                 colour_vector = NULL, feature_colname = "Cell.Type",
                                 cex = 1, layered = FALSE) {
    # if plotting the structure, users do not have to enter the params
    # we have stored the categories and colours for them
    if (feature_colname == "Structure" & is.null(categories_of_interest)) {
        categories_of_interest <- c("Border",
                                    "Inside", 
                                    "Infiltrated.CoI",
                                    "Outside",    
                                    "Stromal.CoI", 
                                    "Internal.margin",     
                                    "Internal.margin.CoI",
                                    "External.margin", 
                                    "External.margin.CoI")
        colour_vector <- c("black", "pink", "purple", "yellow", "orange", "lightgreen", "darkgreen", "lightblue", "blue")
    }
    
    # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
    Cell.X.Position <- Cell.Y.Position <- Category <- NULL
    
    if (methods::is(spe_object, 'SpatialExperiment')){
        # formatted_data <- data.frame(SummarizedExperiment::colData(spe_object))
        formatted_data <- get_colData(spe_object)
    }
    else formatted_data <- spe_object
    
    #CHECK
    if (length(categories_of_interest) != length(colour_vector)) {
        stop("The colour vector is not the same length as the phenotypes of interest")
    }
    
    # if some categories are not in the data, delete them from the categories_of_interest vector
    # delete the corresponding colours as well
    # return a message informing the deleted category
    for (category in categories_of_interest) {
        if (!(category %in% unique(formatted_data[[feature_colname]]))) {
            #stop(paste(category, "cells were not found"), sep="")
            cat_idx <- match(category, categories_of_interest)
            categories_of_interest <- categories_of_interest[-cat_idx]
            colour_vector <- colour_vector[-cat_idx]
            methods::show(paste(category, "cells were not found and not plotted", sep=" "))
        }
    }
    
    #set all categories of those that aren't in categories_of_interest to be "OTHER"
    if (any(!formatted_data[[feature_colname]] %in% categories_of_interest)) {
        formatted_data[!formatted_data[[feature_colname]] %in% categories_of_interest,
        ][[feature_colname]] <- "OTHER"
    }
    
    #Assign the colour to corresponding phenotypes in df
    formatted_data$color <- ""
    for (category in categories_of_interest) {
        idx <- which(categories_of_interest == category)
        formatted_data[formatted_data[[feature_colname]] == category, ]$color <- colour_vector[idx]
    }
    
    if (any(formatted_data[[feature_colname]] == "OTHER")) {
        formatted_data[formatted_data[[feature_colname]] == "OTHER", ]$color <- "lightgrey"
        all_categories <- c(categories_of_interest, "OTHER")
        all_colours <- c(colour_vector, "lightgrey")
    } else {
        all_categories <- categories_of_interest
        all_colours <- colour_vector
    }
    
    if (layered){
        all_cell_types_ordered <- c(categories_of_interest, 
                                    setdiff(unique(formatted_data[[feature_colname]]), categories_of_interest))
        formatted_data[[feature_colname]] <- as.factor(formatted_data[[feature_colname]])
        levels(formatted_data[[feature_colname]]) <- all_cell_types_ordered
        
        p <- ggplot(formatted_data, aes(x = .data$Cell.X.Position, y = .data$Cell.Y.Position)) +
            geom_point(aes(colour = .data[[feature_colname]]), size = cex)
        p <- p +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"))+
            guides(alpha = "none") +
            ggtitle(paste(attr(spe_object, "name"), feature_colname, sep = " ")) +
            scale_color_manual(breaks = categories_of_interest, values=colour_vector)
    }
    else{
        p <- ggplot(formatted_data, aes(x = .data$Cell.X.Position, y = .data$Cell.Y.Position)) +
            geom_point(aes(colour = .data[[feature_colname]]), size = cex)
        p <- p +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"))+
            guides(alpha = "none") +
            ggtitle(paste(attr(spe_object, "name"), feature_colname, sep = " ")) +
            scale_color_manual(breaks = all_categories, values=all_colours)
    }
    return(p)
}

