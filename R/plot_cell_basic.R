#' plot_cell_basic
#' 
#' @description Produces a scatter plot of the cells in the tissue. Cells are coloured
#' categorically by specified column. Cells not part of the celltypes of interest will be coloured "lightgrey"
#' 
#' @param sce_object Singlecellexperiment object in the form of the output of format_image_to_sce
#' @param types_of_interest Vector of cell cell types to be coloured
#' @param colour_vector Vector specifying the colours of each cell phenotype
#' @param column String specifying the column to be coloured
#' @param cex Number specifying the size of the points on the plot. Default is 0.4.
#' @import dplyr
#' @importFrom SummarizedExperiment colData assay
#' @importFrom graphics legend par 
#' @return A plot is returned
#' @export

plot_cell_basic <- function (sce_object, types_of_interest, colour_vector, column, cex = 0.4) {
  
  
  Cell.X.Position <- Cell.Y.Position <- NULL
  assign(column, NULL)
  formatted_data <- data.frame(colData(sce_object))

  if (length(types_of_interest) != length(colour_vector)) {
    stop("The colour vector is not the same length as the celltypes of interest")
  }
  
  real_celltypes <- types_of_interest
  for (phenotype in types_of_interest) {
    if (!(phenotype %in% unique(formatted_data[[column]]))) {
      print(paste(phenotype, "cells were not found"), sep = "")
      real_celltypes <- real_celltypes[real_celltypes != phenotype]
    }
  }
  
  colour_vector <- colour_vector[match(real_celltypes, types_of_interest)]
  
  types_of_interest <- real_celltypes
  
  if (any(!formatted_data[[column]] %in% types_of_interest)) {
    formatted_data[!formatted_data[[column]] %in% types_of_interest, 
    ][[column]] <- "OTHER"
  }
  
  formatted_data$color <- ""
  
  for (phenotype in types_of_interest) {
    idx <- which(types_of_interest == phenotype)
    formatted_data[formatted_data[[column]] == phenotype, 
    ]$color <- colour_vector[idx]
  }
  
  if (any(formatted_data[[column]] == "OTHER")) {
    formatted_data[formatted_data[[column]] == "OTHER", ]$color <- "grey"
    all_phenotypes <- c(types_of_interest, "OTHER")
    all_colours <- c(colour_vector, "grey")
  }
  
  else {
    all_phenotypes <- types_of_interest
    all_colours <- colour_vector
  }

  name_of_object <- attr(sce_object, "name")
  
  plot(formatted_data$Cell.X.Position,formatted_data$Cell.Y.Position,
       pch = 19,cex = cex, col = formatted_data$color, 
       xlab = "X Position", ylab = "Y Position",
       main = paste("Plot", name_of_object, "by" ,column))
  
  x.max <- max(sce_object$Cell.X.Position)
  y.max <- max(sce_object$Cell.Y.Position)
  x <- x.max
  y <- y.max/2
  par(mar=c(5.1, 4.1, 4.1, 5.1), xpd=TRUE)
  legend(x,y, legend=all_phenotypes,
         col=all_colours, cex=0.6, pch = 19,box.lty=0, bg='gray')
  
}
