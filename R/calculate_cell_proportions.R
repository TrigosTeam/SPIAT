#' calculate_cell_proportions
#'
#' @description Calculate the number and proportion of each cell phenotype in an image
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param reference_celltypes Vector specifying reference phenotypes. If NULL (default), then the proportion of each cell type against all cells is returned.
#' Alternatively, a custom vector of cell types can be used as input, and these will be used as the denominator in the calculation of the proportions.
#' @param celltypes_to_exclude Vector specifying cell types to exclude. For example "OTHER" will exclude that celltype from the Total. If NULL, all cell types are included
#' @param column Column of cells to choose the phenotype from (e.g. Phenotype, Cell.Type, etc)
#' @importFrom SummarizedExperiment colData
#' @importFrom stats complete.cases
#' @return A data.frame is returned
#' @examples
#' @export

calculate_cell_proportions <- function(sce_object, column="Phenotype", reference_celltypes = NULL, celltypes_to_exclude = NULL){

    #Reads the image file and deletes cell rows with NA positions
    cell_loc <- data.frame(colData(sce_object))
    cell_loc <- cell_loc[complete.cases(cell_loc),]
    
    #CHECK
    if (nrow(cell_loc) == 0) {
        stop("No cells found for calculating cell proportions")
    }

    #Creates frequency/bar plot of all cell types in the entire image
    cell_proportions <- as.data.frame(table(cell_loc[,column]))
    names(cell_proportions)[1] <- 'Cell_type'
    names(cell_proportions)[2] <- 'Number_of_celltype'
    
    #Exclude any phenotypes not wanted
    if (!is.null(celltypes_to_exclude)) {
        for (celltype in celltypes_to_exclude) {
            cell_proportions <- cell_proportions[!grepl(celltype, cell_proportions$Cell_type), ]
        }
    } 
    
    if(is.null(reference_celltypes)){
        celltype_cells <- cell_proportions
        celltype_cells_total <- sum(celltype_cells$Number_of_celltype)
        celltype_cells$Proportion <- celltype_cells$Number_of_celltype/celltype_cells_total
        celltype_cells$Percentage <- celltype_cells$Proportion*100
        celltype_cells$Proportion_name <- paste(celltype_cells$Cell_type, "Total", sep="/")  
        cell_proportions$Reference <- "Total"
        cell_proportions <- celltype_cells
    
    }else{
        celltype_cells <- cell_proportions
        celltype_cells_total <- sum(celltype_cells$Number_of_celltype[celltype_cells$Cell_type %in% reference_celltypes])
        celltype_cells$Proportion <- celltype_cells$Number_of_celltype/celltype_cells_total
        celltype_cells$Percentage <- celltype_cells$Proportion*100
        celltype_cells$Proportion_name <- paste(celltype_cells$Cell_type, "Custom", sep="/")   
        cell_proportions <- celltype_cells
        cell_proportions$Reference <- paste(reference_celltypes, collapse=",")
    }
    
    #order by Reference celltype (reverse to have Total first if present) then by highest proportion
    cell_proportions <- cell_proportions[rev(order(cell_proportions$Proportion)), ]
    
    return(cell_proportions)
}
