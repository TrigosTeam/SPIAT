#' plot_cell_categories
#' 
#' @description Produces a scatter plot of the cells in the tissue. Cells are coloured
#' categorically by phenotype. Cells not part of the phenotypes of interest will be coloured "lightgrey"
#' 
#' @param sce_object Singlecellexperiment object in the form of the output of format_image_to_sce
#' @param phenotypes_of_interest Vector of cell phenotypes to be coloured
#' @param colour_vector Vector specifying the colours of each cell phenotype
#' @import dplyr
#' @import ggplot2
#' @importFrom SummarizedExperiment colData assay
#' @return A plot is returned
#' @examples
#' phenotypes_of_interest <- c("AMACR", "CD3,CD8", "PDL-1")
#' colour_vector <- c("darkgrey", "blue", "red")
#' plot_cell_categories(SPIAT::formatted_image, phenotypes_of_interest, colour_vector)
#' @export

plot_cell_categories <- function(sce_object, phenotypes_of_interest, colour_vector) {
  
  set.seed(2020)
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Cell.X.Position <- Cell.Y.Position <- Phenotype <- NULL
  
  formatted_data <- data.frame(colData(sce_object))
  
  formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column
  
  intensity_matrix <- assay(sce_object)
  
  markers <- rownames(intensity_matrix)
  cell_ids <- colnames(intensity_matrix)
  
  rownames(intensity_matrix) <- NULL
  colnames(intensity_matrix) <- NULL
  intensity_matrix_t <- t(intensity_matrix)
  intensity_df <- data.frame(intensity_matrix_t)
  colnames(intensity_df) <- markers
  
  formatted_data <- cbind(formatted_data, intensity_df)
  formatted_data <- formatted_data[complete.cases(formatted_data),]
  
  #CHECK
  if (length(phenotypes_of_interest) != length(colour_vector)) {
    stop("The colour vector is not the same length as the phenotypes of interest")
  }
  for (phenotype in phenotypes_of_interest) {
    if (!(phenotype %in% unique(formatted_data$Phenotype))) {
      stop(paste(phenotype, "cells were not found"), sep="")
    }
  }
  
  #set all phenotypes of those that aren't in phenotypes_of_interest to be "OTHER"
  if (any(!formatted_data$Phenotype %in% phenotypes_of_interest)) {
    formatted_data[!formatted_data$Phenotype %in% phenotypes_of_interest,]$Phenotype <- "OTHER"
  }
  
  
  #Assign the colour to corresponding phenotypes in df
  formatted_data$color <- ""
  for (phenotype in phenotypes_of_interest) {
    idx <- which(phenotypes_of_interest == phenotype)
    formatted_data[formatted_data$Phenotype == phenotype, ]$color <- colour_vector[idx]
  }
  if (any(formatted_data$Phenotype == "OTHER")) {
    formatted_data[formatted_data$Phenotype == "OTHER", ]$color <- "lightgrey"
    all_phenotypes <- c(phenotypes_of_interest, "OTHER")
    all_colours <- c(colour_vector, "lightgrey")
  } else {
    all_phenotypes <- phenotypes_of_interest
    all_colours <- colour_vector
  }
  
  p <- ggplot(formatted_data, aes(x = Cell.X.Position, y = Cell.Y.Position, colour = Phenotype))
  if (any(formatted_data$Phenotype == "OTHER")) {
    p <- p + geom_point(data=subset(formatted_data, Phenotype=='OTHER'), aes(colour = Phenotype), size = 1) + 
      geom_point(data=subset(formatted_data, Phenotype!='OTHER'), aes(colour = Phenotype), size = 1) 
  }
  else {
    p <- p + geom_point(aes(colour = Phenotype), size = 1)
  }
  p <- p +
    guides(alpha = FALSE) +
    labs(colour = "Phenotypes") + 
    scale_color_manual(breaks = all_phenotypes, values=all_colours) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  print(p)
}
