#' calculate_cell_proportions
#'
#' @description Calculates the number and proportion of each cell type.
#'
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param reference_celltypes String Vector specifying reference cell types. If
#'   NULL (default), then the proportion of each cell type against all cells is
#'   returned. Alternatively, a custom vector of cell types can be used as
#'   input, and these will be used as the denominator in the calculation of the
#'   proportions.
#' @param celltypes_to_exclude String Vector specifying cell types to exclude.
#'   For example "OTHER" will exclude that celltype from the Total. If NULL, all
#'   cell types are included.
#' @param feature_colname String. Column of cells to choose the cell type from
#'   (e.g. Phenotype, Cell.Type, etc).
#' @param plot.image Boolean. Whether to plot the barplot of the cell percentages.
#'   By default is TRUE.
#' @import ggplot2
#' @return A data.frame is returned
#' @examples
#' calculate_cell_proportions(SPIAT::defined_image, reference_celltypes = NULL, 
#' celltypes_to_exclude = "Others", feature_colname="Cell.Type", plot.image = FALSE)
#' @export

calculate_cell_proportions <- function(spe_object,  reference_celltypes = NULL, 
                                       celltypes_to_exclude = NULL, 
                                       feature_colname="Phenotype",plot.image = TRUE){
    Cell_type <- Percentage <- NULL
    #Reads the image file and deletes cell rows with NA positions
    cell_loc <- get_colData(spe_object)
    
    #CHECK
    if (nrow(cell_loc) == 0) {
        stop("No cells found for calculating cell proportions")}

    #Creates frequency/bar plot of all cell types in the entire image
    cell_proportions <- as.data.frame(table(cell_loc[,feature_colname]))
    names(cell_proportions)[1] <- 'Cell_type'
    names(cell_proportions)[2] <- 'Number_of_celltype'
    
    #Exclude any phenotypes not wanted
    if (!is.null(celltypes_to_exclude)) {
        for (celltype in celltypes_to_exclude) {
            cell_proportions <- 
                cell_proportions[!grepl(celltype, 
                                        cell_proportions[["Cell_type"]]), ]
        }} 
    if(is.null(reference_celltypes)){
        celltype_cells <- cell_proportions
        celltype_cells_total <- sum(celltype_cells$Number_of_celltype)
        celltype_cells$Proportion <- celltype_cells$Number_of_celltype/celltype_cells_total
        celltype_cells$Percentage <- celltype_cells$Proportion*100
        celltype_cells$Proportion_name <- paste(celltype_cells[[feature_colname]],
                                                "Total", sep="/")  
        cell_proportions$Reference <- "Total"
        cell_proportions <- celltype_cells
    
    }else{
        celltype_cells <- cell_proportions
        celltype_cells_total <- 
            sum(celltype_cells$Number_of_celltype[celltype_cells[[feature_colname]] 
                                                  %in% reference_celltypes])
        celltype_cells$Proportion <- celltype_cells$Number_of_celltype/celltype_cells_total
        celltype_cells$Percentage <- celltype_cells$Proportion*100
        celltype_cells$Proportion_name <- paste(celltype_cells[[feature_colname]], 
                                                "Custom", sep="/")   
        cell_proportions <- celltype_cells
        cell_proportions$Reference <- paste(reference_celltypes, collapse=",")
    }
    
    #order by Reference celltype (reverse to have Total first if present) then by highest proportion
    cell_proportions <- cell_proportions[rev(order(cell_proportions$Proportion)), ]
    
    if(plot.image){
        g <- ggplot(cell_proportions, aes(x=Cell_type, y=Percentage)) +
            geom_bar(stat='identity') + theme_bw()
        methods::show(g)}
    
    return(cell_proportions)
}
