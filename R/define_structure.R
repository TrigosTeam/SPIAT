#' define_structure
#'
#' @description After identifying the bordering cells of tissue regions and
#'   calculating the distances of each cell to the bordering cells, this
#'   function further identifies the cells that are located in the inside and
#'   outside of the identified regions, and in the internal and external
#'   margins. It also identifies particular types of cells that are infiltrated,
#'   stromal, internal margin or external margin cells.
#'
#' @param spe_object SpatialExperiment object that contains information of
#'   tumour bordering cells and cell distances to border (`colData()` has
#'   `Region` and `Distance.To.Border` columns).
#' @param cell_types_of_interest String Vector of the names of the particular
#'   types of cells.
#' @param n_margin_layers Integer. The number of layers of cells that compose
#'   the internal/external margins.
#' @param feature_colname String Specifying the column that contains the names
#'   of the immune cells.
#' @import dplyr
#' @export
#' @return A new spe object is returned. Under the `Region` column, there will
#'   be potential categories including `Border` - the bordering cells,
#'   `Infiltrated.CoI` - cells of interest that present inside of the tissue
#'   regions, `Inside` - cells within the regiona excluding the
#'   `Infiltrated.CoI` cells and the cells at internal margin, `Stromal.CoI` -
#'   cells of interest that present outside of the tissue regions, `Outside` -
#'   cells outside of the tissue regions excluding the `Stromal.CoI` cells,
#'   `Internal.margin.CoI` - cells of interest that are in the internal margin
#'   of the tissue regions, `Internal.margin` - cells in the internal margin of
#'   the tissue regions excluding the `Internal.margin.CoI` cells,
#'   `External.margin.CoI` - cells of interest that are in the external margin
#'   of the tissue regions, `External.margin` - cells in the external margin of
#'   the tissue regions excluding the `External.margin.CoI` cells.
#' @examples
#' spe_border <- identify_bordering_cells(SPIAT::defined_image,
#' reference_cell = "Tumour", feature_colname = "Cell.Type", n_to_exclude = 10)
#' spe_dist <- calculate_distance_to_margin(spe_border)
#' spe_structure <- define_structure(spe_dist,
#' cell_types_of_interest = c("Immune1","Immune2","Immune3"),
#' feature_colname = "Cell.Type", n_margin_layers = 5)
#' plot_cell_categories(spe_structure, feature_colname = "Structure")

define_structure <- function(spe_object, cell_types_of_interest, 
                             feature_colname = "Cell.Type",
                             n_margin_layers = 5){
    
    # calculate the width of internal/external margin
    min_dist <- average_minimum_distance(spe_object)
    margin_dist <- n_margin_layers * min_dist
    sprintf("How many layers of cells in the external/internal margin: %i", 
            n_margin_layers)
    sprintf("The width of internal/external margin: %f", margin_dist)
    
    #CHECK if the distance to bordering cells is calculated
    if (is.null(spe_object$Distance.To.Border)){
        stop(sprintf("Distance.To.Border not calculated yet for %i", 
                     deparse(substitute(spe_object))))
    }
    
    data <- data.frame(SummarizedExperiment::colData(spe_object))
    data[,"Structure"] <- data$Region
    data[intersect(which(data$Region == "Inside"),
                   which(data[[feature_colname]] %in% cell_types_of_interest)), 
         "Structure"] <- "Infiltrated.CoI"
    data[intersect(which(data$Region == "Outside"),
                   which(data[[feature_colname]] %in% cell_types_of_interest)), 
         "Structure"] <- "Stromal.CoI"
    data[intersect(which(data$Distance.To.Border < margin_dist), 
                   which(data$Region == "Inside")), 
         "Structure"] <- "Internal.margin"
    data[intersect(which(data$Structure == "Internal.margin"), 
                   which(data[[feature_colname]] %in% cell_types_of_interest)), 
         "Structure"] <- "Internal.margin.CoI"
    data[intersect(which(data$Distance.To.Border < margin_dist), 
                   which(data$Region == "Outside")), 
         "Structure"] <- "External.margin"
    data[intersect(which(data$Structure == "External.margin"), 
                   which(data[[feature_colname]] %in% cell_types_of_interest)), 
         "Structure"] <- "External.margin.CoI"
    
    SummarizedExperiment::colData(spe_object)$Structure <- data[,"Structure"]
    
    return(spe_object)
}
