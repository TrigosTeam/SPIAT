#' measure_association_to_cell_properties
#'
#' @description Plots the density or boxplot of a property of two cell phenotypes or compares using t test/wilcoxon rank sum test
#' 
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param property String that is the name of the column of interest 
#' @param phenotypes Vector of phenotypes of interest
#' @param merge Vector of phenotypes to be merged
#' @param merge_name String that is the name of the merged phenotype
#' @param method the analysis to do on the selected phenotypes and property. Options are density, box, t, wilcox 
#' @param Nucleus.Ratio when the ratio of the nucleus size is of interest
#' @param log.scale if log the data
#' @importFrom SummarizedExperiment colData
#' @importFrom stats t.test wilcox.test
#' @importFrom dittoSeq dittoColors  
#' @import ggplot2
#' @export

measure_association_to_cell_properties <- function(sce_object, property = "Cell.Area", 
                                           phenotypes, merge = NULL, merge_name = NULL, 
                                           method = "density", Nucleus.Ratio = FALSE,
                                           log.scale = FALSE) {

  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Phenotype <- NULL
  
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
    formatted_data["Nucleus.Ratio"] <- formatted_data$Nucleus.Area/formatted_data$Cell.Area
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
    
    # get colourblind-friendly colours
    colours <- dittoColors()[1:2]
    
    p <- ggplot(formatted_data, aes(x=formatted_data[,property], color = Phenotype)) + 
      geom_density() + 
      labs(x = property) +
      scale_color_manual(values = colours) +
      theme_bw()
  }
  
  # Plot the boxplot 
  if (method == "box"){
    p <- ggplot(formatted_data, aes(Phenotype,formatted_data[,property])) + 
      geom_boxplot()  +
      stat_boxplot(coef=3) +
      ylab(property) +
      theme_bw()

  }
  
  if (method == "t"){
    # CHECK
    if (length(phenotypes) != 2){
      stop("wrong number of inputs to do t.test")
    }
    else{
      p <- t.test(formatted_data[formatted_data$Phenotype == phenotypes[1],property],
                  formatted_data[formatted_data$Phenotype == phenotypes[2],property])
      
      # make output name nicer
      p$data.name <- paste(phenotypes[1], "and", phenotypes[2])
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
