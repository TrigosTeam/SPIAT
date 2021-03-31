#' plot_cell_basic
#' 
#' @description Produces a scatter plot of the cells in the tissue. Cells are coloured
#' categorically by Ceell.Type. Cells not part of the celltypes of interest will be coloured "lightgrey"
#' 
#' @param sce_object Singlecellexperiment object in the form of the output of format_image_to_sce
#' @param celltypes_of_interest Vector of cell cell types to be coloured
#' @param colour_vector Vector specifying the colours of each cell phenotype
#' @param column String specifying the column to be coloured
#' @import dplyr
#' @importFrom SummarizedExperiment colData assay
#' @return A plot is returned
#' @examples
#' celltypes_of_interest <- c("AMACR", "CD3,CD8", "PDL-1")
#' colour_vector <- c("darkgrey", "blue", "red")
#' plot_cell_basic(SPIAT::formatted_image, celltypes_of_interest, colour_vector, column)
#' @export

plot_cell_basic <- function (sce_object, celltypes_of_interest, colour_vector, column) 
{
  Cell.X.Position <- Cell.Y.Position <- NULL
  assign(column, NULL)
  formatted_data <- data.frame(colData(sce_object))

  
  if (length(celltypes_of_interest) != length(colour_vector)) {
    stop("The colour vector is not the same length as the celltypes of interest")
  }
  
  for (phenotype in celltypes_of_interest) {
    if (!(phenotype %in% unique(formatted_data[[column]]))) {
      stop(paste(phenotype, "cells were not found"), sep = "")
    }
  }
  if (any(!formatted_data[[column]] %in% celltypes_of_interest)) {
    formatted_data[!formatted_data[[column]] %in% celltypes_of_interest, 
    ][[column]] <- "OTHER"
  }
  formatted_data$color <- ""
  for (phenotype in celltypes_of_interest) {
    idx <- which(celltypes_of_interest == phenotype)
    formatted_data[formatted_data[[column]] == phenotype, 
    ]$color <- colour_vector[idx]
  }
  if (any(formatted_data[[column]] == "OTHER")) {
    formatted_data[formatted_data[[column]] == "OTHER", ]$color <- "grey"
    all_phenotypes <- c(celltypes_of_interest, "OTHER")
    all_colours <- c(colour_vector, "grey")
  }
  else {
    all_phenotypes <- celltypes_of_interest
    all_colours <- colour_vector
  }
  
  plot(formatted_data$Cell.X.Position,formatted_data$Cell.Y.Position,
            pch = 19,cex = 0.2, col = formatted_data$color)
  x.max <- max(sce_object$Cell.X.Position)
  y.max <- max(sce_object$Cell.Y.Position)
  x <- x.max
  y <- y.max/2
  par(mar=c(5.1, 4.1, 4.1, 5.1), xpd=TRUE)
  legend(x,y, legend=all_phenotypes,
         col=all_colours, cex=0.6, pch = 19,box.lty=0, bg='gray')
  
}