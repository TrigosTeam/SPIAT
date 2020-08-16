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
#' @export

plot_cell_categories <- function(sce_object, phenotypes_of_interest, colour_vector) {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Cell.X.Position <- Cell.Y.Position <- Phenotype <- NULL
  
  formatted_data <- data.frame(colData(sce_object))
  
  formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column
  
  expression_matrix <- assay(sce_object)
  
  markers <- rownames(expression_matrix)
  cell_ids <- colnames(expression_matrix)
  
  rownames(expression_matrix) <- NULL
  colnames(expression_matrix) <- NULL
  expression_matrix_t <- t(expression_matrix)
  expression_df <- data.frame(expression_matrix_t)
  colnames(expression_df) <- markers
  
  formatted_data <- cbind(formatted_data, expression_df)
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
  formatted_data[!formatted_data$Phenotype %in% phenotypes_of_interest,]$Phenotype <- "OTHER"
  
  
  #Assign the colour to corresponding phenotypes in df
  formatted_data$color <- ""
  for (phenotype in phenotypes_of_interest) {
    idx <- which(phenotypes_of_interest == phenotype)
    formatted_data[formatted_data$Phenotype == phenotype, ]$color <- colour_vector[idx]
  }
  formatted_data[formatted_data$Phenotype == "OTHER", ]$color <- "lightgrey"
  
  
  all_phenotypes <- c(phenotypes_of_interest, "OTHER")
  all_colours <- c(colour_vector, "lightgrey")
  
  p <- ggplot(formatted_data, aes(x = Cell.X.Position, y = Cell.Y.Position, colour = Phenotype)) +
    geom_point(aes(colour = Phenotype), size = 1) +
    guides(alpha = F) +
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