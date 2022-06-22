#' calculate_proportions_of_cells_in_structure
#'
#' @description Calculate the proportion of cells of interest in each defined
#'   tumour structure relative to all cells in each structure and relative to
#'   the same cell type in the whole image.
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param cell_types_of_interest String Vector of immune cells to consider.
#' @param feature_colname String. The name of the column where the immune cell
#'   types are under.
#' @return A data.frame
#' @export
#' @examples
#' spe_border <- identify_bordering_cells(SPIAT::defined_image,
#' reference_cell = "Tumour", feature_colname = "Cell.Type", n_to_exclude = 10)
#' spe_dist <- calculate_distance_to_tumour_margin(spe_border)
#' spe_structure <- define_structure(spe_dist,
#' names_of_immune_cells = c("Immune1","Immune2","Immune3"),
#' feature_colname = "Cell.Type", n_margin_layers = 5)
#' calculate_proportions_of_cells_in_structure(spe_structure,
#' cell_types_of_interest = c("Immune1","Immune3"),feature_colname="Cell.Type")

calculate_proportions_of_cells_in_structure <- function(spe_object, 
                                                        cell_types_of_interest, 
                                                        feature_colname){
    data_local <- get_colData(spe_object)
    
    #Relative to all cells in each area
    proportions <- vector()
    for(cell_type in cell_types_of_interest){
        temp <- data_local[data_local[,feature_colname] == cell_type,]
        p_i_in <- length(which(temp$Structure == "Infiltrated.immune"))/
            length(which(data_local$Structure == "Inside"))
        p_i_i.f.in <- length(which(temp$Structure =="Internal.margin.immune"))/
            length(which(data_local$Structure == "Internal.margin"))
        p_i_i.f.out <- length(which(temp$Structure=="External.margin.immune"))/
            length(which(data_local$Structure == "External.margin"))
        p_i_out <- length(which(temp$Structure == "Stromal.immune"))/
            length(which(data_local$Structure == "Outside"))
        
        proportions <- rbind(proportions,
                             c(cell_type, "All_cells_in_the_structure", 
                               p_i_in, p_i_i.f.in, p_i_i.f.out, p_i_out))
    }
    
    proportions <- as.data.frame(proportions)
    
    #Relative to all interested cells in each area
    data_local_interested <- data_local[data_local[,feature_colname] %in% 
                                            cell_types_of_interest,]
    for(cell_type in cell_types_of_interest){
        temp <- data_local[data_local[,feature_colname] == cell_type,]
        p_i_in <- length(which(temp$Structure == "Infiltrated.immune"))/
            length(which(data_local_interested$Structure == 
                             "Infiltrated.immune"))
        p_i_i.f.in <- length(which(temp$Structure == "Internal.margin.immune"))/
            length(which(data_local_interested$Structure == 
                             "Internal.margin.immune"))
        p_i_i.f.out <- length(which(temp$Structure=="External.margin.immune"))/
            length(which(data_local_interested$Structure == 
                             "External.margin.immune"))
        p_i_out <- length(which(temp$Structure == "Stromal.immune"))/
            length(which(data_local_interested$Structure == "Stromal.immune"))
        
        proportions <- rbind(
            proportions,c(cell_type, "All_cells_of_interest_in_the_structure", 
                          p_i_in, p_i_i.f.in, p_i_i.f.out, p_i_out))
    }
    
    #Relative to same type of cell in the whole image
    for(cell_type in cell_types_of_interest){
        temp <- data_local[data_local[,feature_colname] == cell_type,]
        p_i_in <- length(which(temp$Structure == "Infiltrated.immune"))/
            dim(temp)[1]
        p_i_i.f.in <- length(which(temp$Structure == "Internal.margin.immune"))/
            dim(temp)[1]
        p_i_i.f.out <- length(which(temp$Structure =="External.margin.immune"))/
            dim(temp)[1]
        p_i_out <- length(which(temp$Structure == "Stromal.immune"))/
            dim(temp)[1]
        
        proportions <- rbind(
            proportions, c(cell_type, "The_same_cell_type_in_the_whole_image", 
                           p_i_in, p_i_i.f.in, p_i_i.f.out, p_i_out))
    }
    
    #Total interested cells relative to all cells in each area
    p_i_in <- length(which(data_local$Structure == "Infiltrated.immune"))/
        length(which(data_local$Structure == "Inside"))
    p_i_i.f.in <- length(which(data_local$Structure=="Internal.margin.immune"))/
        length(which(data_local$Structure == "Internal.margin"))
    # proportion of immune cells of the invasive front that are outside
    p_i_i.f.out <-length(which(data_local$Structure=="External.margin.immune"))/
        length(which(data_local$Structure == "External.margin"))
    # proportion of exclusice immune cells that are out of the tumour region
    p_i_out <- length(which(data_local$Structure == "Stromal.immune"))/
        length(which(data_local$Structure == "Outside"))
    
    proportions <- rbind(
        proportions, c("All_cells_of_interest", "All_cells_in_the_structure", 
                       p_i_in, p_i_i.f.in, p_i_i.f.out, p_i_out))
    
    
    colnames(proportions) <- c("Cell.Type", "Relative_to", 
                               "P.Infiltrated.Immune",
                               "P.Internal.Margin.Immune",
                               "P.External.Margin.Immune", "P.Stromal.Immune")
    proportions$P.Infiltrated.Immune <- 
        as.numeric(proportions$P.Infiltrated.Immune)
    proportions$P.Stromal.Immune <- as.numeric(proportions$P.Stromal.Immune)
    proportions$P.Internal.Margin.Immune <- 
        as.numeric(proportions$P.Internal.Margin.Immune)
    proportions$P.External.Margin.Immune <- 
        as.numeric(proportions$P.External.Margin.Immune)
    
    return(proportions)
}
