#' calculate_cell_proportions
#'
#' @description Calculate the number and proportion of each cell phenotype in image
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param reference_celltypes Vector specifying reference phenotypes. For example, "CD3" will calculate the proportion
#' of each CD3-containing phenotype against all CD3. "Total" can be used to calculate the proportion of each phenotype against total cells and is the default.
#' @param celltypes_to_exclude Vector specifying celltypes to exclude. For example "OTHER" will exclude that celltype from the Total (optional)
#' @param column Column of cells to choose the phenotype from (e.g. Cell.Type, Cell.Type2, etc)
#' @importFrom SummarizedExperiment colData
#' @importFrom stats complete.cases
#' @return A data.frame is returned
#' @examples
#' p_cells <- calculate_cell_proportions(SPIAT::formatted_image, reference_celltypes=c("Total", "CD3"))
#' @export

calculate_cell_proportions <- function(sce_object, column="Phenotype", reference_celltypes = "Total", celltypes_to_exclude = NULL){

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
    
    #Calculate proportions for each celltype specified   
    datalist = list()
    
    for (celltype in reference_celltypes) {
        
        if (celltype != "Total") {
            
            celltype_cells <- cell_proportions[grepl(celltype, cell_proportions$Cell_type), ]
            
        } else {
        
            celltype_cells <- cell_proportions
            
        }
        
            celltype_cells_total <- sum(celltype_cells$Number_of_celltype)
            celltype_cells$Reference <- rep(celltype, nrow(celltype_cells))
            celltype_cells$Number_of_reference <- rep(celltype_cells_total, nrow(celltype_cells))
            celltype_cells$Proportion <- celltype_cells$Number_of_celltype/celltype_cells_total
            celltype_cells$Percentage <- celltype_cells$Proportion*100
            celltype_cells$Proportion_name <- paste(celltype_cells$Cell_type, celltype, sep="/")
        
        datalist[[celltype]] <- celltype_cells
        
    }
    
    cell_proportions <- do.call(rbind, datalist)
    rownames(cell_proportions) <- NULL
    
    #order by Reference celltype (reverse to have Total first if present) then by highest proportion
    cell_proportions <- cell_proportions[rev(order(cell_proportions$Reference, cell_proportions$Proportion)), ]
    
    
    return(cell_proportions)
}
