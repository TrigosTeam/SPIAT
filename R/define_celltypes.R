#' define_celltypes
#'
#' @description Define new cell types based on the existing cell types
#'   (categories) under a selected column (e.g. base on marker combinations
#'   under "Phenotype" column). This function will create a new column to store
#'   the new cell types.
#'
#' @details Users need to specify the names of the old cell categories and under
#'   which column the old cell categories exist. Then the users specify the
#'   names of the new cell types and the name of the new column to store the new
#'   cell types. Any cell categories that are not specified in `categories` arg
#'   but present in the image will be defined as "Undefined" in the new column.
#'
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param categories Vector. Names of the old cell types to be defined; if NULL,
#'   the function will use predefined categories and names
#' @param category_colname (Phenotype) String specifying the name of the column
#'   having the categories to be defined, by default "Phenotype".
#' @param names Vector of new names assigned to the selected categories; if
#'   NULL, the function will use predefined categories and names. Should be of
#'   the same length of `categories`.
#' @param new_colname (Optional) String specifying the name of the column to be
#'   added, by default "Cell.Type".
#' @param print_names (Optional) Boolean if the user wants the original and new
#'   names printed. Default is FALSE.
#' @export
#' @return An new SPE object is returned
#' @examples
#' # the selected column is:
#' category_colname = "Phenotype"
#' # define the following marker combinations:
#' categories <- c("Tumour_marker", "Immune_marker1,Immune_marker2", 
#' "Immune_marker1,Immune_marker3",
#' "Immune_marker1,Immune_marker2,Immune_marker4", "OTHER")
#' # the new defined cell names:
#' names = c("Tumour", "Immune1", "Immune2","Immune3", "Others")
#' # the new names are stored under this column:
#' new_colname <- "Cell.Type"
#'
#' defined_spe <- define_celltypes(SPIAT::simulated_image, 
#' categories = categories, category_colname = category_colname, names = names, 
#' new_colname = new_colname)

define_celltypes <- function(spe_object,categories = NULL, 
                             category_colname = "Phenotype", names = NULL, 
                             new_colname = "Cell.Type", print_names = FALSE){
    
    # CHECK
    if (length(categories) != length(names)){
        methods::show("`length(categories) != length(names)`")
        stop("The old and new cell type names should be of the same length!")
    }
    # default setting
    if (is.null(categories) & is.null(names)){
        pre_categories <- c("AMACR", "CD3,CD4", "CD3,CD8", "CD3", 
                            "AMACR,PDL-1","PDL-1")
        pre_names <- c("Tumour","CD4","CD8","CD3","Tumour,PDL1","PDL1")
    }
    
    else{
        pre_categories <- categories
        pre_names <- names
    }
    
    all_categories <- 
        unique(SummarizedExperiment::colData(spe_object)[[category_colname]])
    categories <- intersect(all_categories, pre_categories)
    names <- pre_names[match(categories, pre_categories)]
    
    if (print_names){
        methods::show(paste("Define new cell types basing on:", categories))
        methods::show(paste("The new cell types are:",names))
    }
    
    spe_object[[new_colname]] <- ""
    names<- c(names,rep("Undefined",length(all_categories[
        !(all_categories %in% categories)])))
    categories <- c(categories,all_categories[
        !(all_categories %in% categories)])
    for (i in seq_len(length(categories))){
        spe_object[,spe_object[[category_colname]] == 
                       categories[i]][[new_colname]] <- names[i]
    }
    
    return(spe_object)
}

