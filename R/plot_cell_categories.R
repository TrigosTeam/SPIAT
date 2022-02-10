#' plot_cell_categories
#' 
#' @description Produces a scatter plot of the cells in the tissue. Cells are coloured
#' categorically by phenotype. Cells not part of the phenotypes of interest will be coloured "lightgrey"
#' 
#' @param sce_object Singlecellexperiment object in the form of the output of format_image_to_sce
#' @param categories_of_interest Vector of cell categories to be coloured
#' @param colour_vector Vector specifying the colours of each cell phenotype
#' @param feature_colname String specifying the column the cell categories belong to
#' @import dplyr
#' @import ggplot2
#' @import tibble
#' @importFrom SummarizedExperiment colData assay
#' @return A plot is returned
#' @examples
#' categories_of_interest <- c("AMACR", "CD3,CD8", "PDL-1")
#' colour_vector <- c("red", "blue", "orange")
#' plot_cell_categories(SPIAT::formatted_image, categories_of_interest, colour_vector)
#' @export

plot_cell_categories <- function(sce_object, categories_of_interest = NULL, 
                                 colour_vector = NULL, feature_colname = "Cell.Type") {
  
  # if plotting the structure, users do not have to enter the params
  # we have stored the categories and colours for them
  if (feature_colname == "Structure" & is.null(categories_of_interest)) {
    categories_of_interest <- c("Border",
                                "Inside", 
                                "Infiltrated.immune",
                                "Outside",    
                                "Stromal.immune", 
                                "Internal.margin",     
                                "Internal.margin.immune",
                                "External.margin", 
                                "External.margin.immune")
    colour_vector <- c("black", "pink", "purple", "yellow", "orange", "lightgreen", "darkgreen", "lightblue", "blue")
  }
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Cell.X.Position <- Cell.Y.Position <- Category <- NULL
  formatted_data <- data.frame(colData(sce_object))
  
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
      print(paste(category, "cells were not found and not plotted"), sep="")
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
  
  p <- ggplot(formatted_data, aes_string(x = "Cell.X.Position", y = "Cell.Y.Position")) +
    geom_point(aes_string(colour = feature_colname), size = 1)
  # p <- ggplot(formatted_data, aes_string(x = "Cell.X.Position", y = "Cell.Y.Position", colour = "color"))
  # if (any(formatted_data[[feature_colname]] == "OTHER")) {
  #   p <- p + geom_point(data=subset(formatted_data, get(feature_colname) =='OTHER'),
  #                       aes_string(colour = "color"), size = 1) +
  #     geom_point(data=subset(formatted_data, get(feature_colname) !='OTHER'),
  #                aes_string(colour = "color"), size = 1)
  # }else{
  #   p <- p + geom_point(aes_string(colour = "color"), size = 1)}
  # 
  p <- p +
    guides(alpha = "none") +
    ggtitle(paste("Plot", attr(sce_object, "name"), feature_colname, sep = " ")) +
    scale_color_manual(breaks = all_categories, values=all_colours)
  # labs(colour = all_categories) +
  
  #   theme(panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.background = element_rect(fill = "white"),
  #         axis.title.x = element_blank(),
  #         axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         axis.title.y = element_blank(),
  #         axis.text.y = element_blank(),
  #         axis.ticks.y = element_blank())
  
  print(p)
}
