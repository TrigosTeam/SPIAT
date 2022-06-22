#' measure_association_to_cell_properties
#'
#' @description Plots the density or boxplot of a property of two cell celltypes
#'   or compares using t test/wilcoxon rank sum test.
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param property String that is the name of the column of interest.
#' @param celltypes String Vector of celltypes of interest.
#' @param feature_colname String that speficies the column of the cell types.
#' @param method String. The analysis to perform on the selected cell types and
#'   property. Options are "density", "box", "t", "wilcox".
#' @param Nucleus.Ratio Boolean whether the ratio of the nucleus size is of
#'   interest.
#' @param log.scale Boolean whether to log the data.
#' @import ggplot2
#' @return With method "box" or "density a plot is returned. With method "t" or
#'   "wilcox", the text output from the test are returned.
#' @examples
#' measure_association_to_cell_properties(image_no_markers,
#'                                       celltypes = c("Tumour", "Immune1"),
#'                                       feature_colname = "Cell.Type",
#'                                       property = "Cell.Size",
#'                                       method = "box")
#' measure_association_to_cell_properties(image_no_markers,
#'                                       celltypes = c("Tumour", "Immune2"),
#'                                       feature_colname="Cell.Type",
#'                                       property = "Cell.Size",
#'                                       method = "t")
#' @export

measure_association_to_cell_properties <- function(spe_object, property = "Cell.Area", 
                                                   celltypes, feature_colname = "Cell.Type",
                                                   method = "density", Nucleus.Ratio = FALSE,
                                                   log.scale = FALSE) {
    
    formatted_data <- get_colData(spe_object)
    
    #CHECK
    if (is.element(property, colnames(formatted_data)) == FALSE) {
        stop("Property of interest not found")
    }
    
    if (!all(celltypes %in% formatted_data[[feature_colname]])) {
        stop("Cell type not found")
    }
    
    # CHECK if nucleus.ratio is the property of interest
    if (Nucleus.Ratio == TRUE){
        formatted_data["Nucleus.Ratio"] <- formatted_data$Nucleus.Area/formatted_data$Cell.Area
        property <- "Nucleus.Ratio"
    }
    
    #Extract interested property and celltypes
    formatted_data <- formatted_data[which(formatted_data[[feature_colname]] %in% celltypes),
                                     c(feature_colname,property)]
    
    
    # CHECK if log the scale
    if (log.scale == TRUE){
        formatted_data[,property] <- log(formatted_data[,property])
    }
    
    # Plot the density 
    if (method == "density"){
        
        # get colourblind-friendly colours
        colours <- dittoSeq::dittoColors()[seq_len(2)]
        
        p <- ggplot(formatted_data, aes_string(x=property, color = feature_colname)) + 
            geom_density() + 
            labs(x = property) +
            scale_color_manual(values = colours) +
            theme_bw()
    }
    
    # Plot the boxplot 
    if (method == "box"){
        
        #code from ggplot2 to show the number in each group
        give.n <- function(x){
            return(c(y = max(x)+1, label = length(x)))
        }
        
        p <- ggplot(formatted_data, aes_string(feature_colname,property, fill = feature_colname)) + 
            geom_boxplot()  +
            stat_summary(fun.data = give.n, geom = "text", vjust = -0.5) +
            ylab(property) +
            theme_bw()
        
    }
    
    if (method == "t"){
        # CHECK
        if (length(celltypes) != 2){
            stop("wrong number of inputs to do t.test")
        }
        else{
            p <- stats::t.test(formatted_data[formatted_data[[feature_colname]] == celltypes[1],property],
                               formatted_data[formatted_data[[feature_colname]] == celltypes[2],property])
            
            # make output name nicer
            p$data.name <- paste(celltypes[1], "and", celltypes[2])
        }
    }
    if (method == "wilcox"){
        # CHECK
        if (length(celltypes) != 2){
            stop("wrong number of inputs to do Wilcoxon Rank Sum test")
        }
        else{
            p <- stats::wilcox.test(formatted_data[formatted_data[[feature_colname]] == celltypes[1],property],
                                    formatted_data[formatted_data[[feature_colname]] == celltypes[2],property])
        }
    }
    
    return (p)
}
