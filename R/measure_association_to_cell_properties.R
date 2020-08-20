#' measure_association_to_cell_properties
#'
#' @description plot the density or boxplot of properties of certain cell phenotypes
#' compare density plots or boxplots of different cell phenotypes
#' t test/wilcoxon rank sum test of a property of two cell phenotypes
#' 
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param property String that is the name of the column of interest 
#' @param phenotypes Vector of phenotypes of interest
#' @param merge Vector of phenotypes to be merged
#' @param merge_name String that is the name of the merged phenotype
#' @param method the analysis to do on the selected phenotypes and property
#' @param Nucleus.Ratio when the ratio of the nucleus size is of interest
#' @param log.scale if log the data
#' @importFrom SummarizedExperiment colData
#' @import ggplot2
#' @export

measure_association_to_cell_properties <- function(sce_object, property = "Cell.Size", 
                                           phenotypes, merge = NULL, merge_name = NULL, 
                                           method = "density", Nucleus.Ratio = FALSE,
                                           log.scale = FALSE) {
  
  formatted_data <- data.frame(colData(sce_object))

  #CHECK
  if (is.element(property, colnames(formatted_data)) == FALSE) {
    stop("Property of interest not found")
  }
  
  if (!all(phenotypes %in% formatted_data$Phenotype)) {
    stop("phenotype not found")
  }
  
  # CHECK if nucleus.ratio is the property of interest
  if (Nucleus.Ratio == TRUE){
    formatted_data["Nucleus.Ratio"] <- formatted_data$Nucleus.Size/formatted_data$Cell.Size
    property <- "Nucleus.Ratio"
  }
  
  #Extract interested property and phenotypes
  formatted_data <- formatted_data[which(formatted_data$Phenotype %in% phenotypes),
                                   c("Phenotype",property)]
  
  # CHECK if there are any columns to be merged
  if (is.null(merge) == FALSE){
    formatted_data[which(formatted_data$Phenotype %in% merge),]$Phenotype <- merge_name
    # update phenotypes of interest
    phenotypes <- c(setdiff(phenotypes, merge), merge_name)
  }
  
  # CHECK if log the scale
  if (log.scale == TRUE){
    formatted_data[,property] <- log(formatted_data[,property])
  }
  
  # Plot the density 
  if (method == "density"){
    p <- ggplot(formatted_data, aes(x=formatted_data[,property], color = Phenotype)) + 
      geom_density() + 
      labs(x = property)
  }
  
  # Plot the boxplot 
  if (method == "box"){
    p <- ggplot(formatted_data, aes(Phenotype,formatted_data[,property], fill = Phenotype)) + 
      geom_boxplot()  +
      stat_boxplot(coef=3) +
      ylab(property)

  }
  
  if (method == "t"){
    # CHECK
    if (length(phenotypes) != 2){
      stop("wrong number of inputs to do t.test")
    }
    else{
      p <- t.test(formatted_data[formatted_data$Phenotype == phenotypes[1],property],
                  formatted_data[formatted_data$Phenotype == phenotypes[2],property])
    }
  }
  if (method == "wilcox"){
    # CHECK
    if (length(phenotypes) != 2){
      stop("wrong number of inputs to do Wilcoxon Rank Sum test")
    }
    else{
      p <- wilcox.test(formatted_data[formatted_data$Phenotype == phenotypes[1],property],
                       formatted_data[formatted_data$Phenotype == phenotypes[2],property])
    }
  }

  return (p)
}
