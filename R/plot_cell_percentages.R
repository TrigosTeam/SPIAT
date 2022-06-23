#' plot_cell_percentages
#'
#' @description Plots cells proportions as barplots.
#' @param cell_proportions Data Frame. Output from
#'   \code{\link{calculate_cell_proportions}}.
#' @param cells_to_exclude String Vector. Markers to exclude.
#' @param cellprop_colname String. Column to use for y axis names. Default is
#'   "Proportion_name".
#' @import ggplot2
#' @import dplyr
#' @return A plot is returned
#' @examples
#' p_cells <- calculate_cell_proportions(SPIAT::simulated_image)
#' plot_cell_percentages(p_cells)
#' @export
plot_cell_percentages <- function(cell_proportions, cells_to_exclude =NULL, 
                                  cellprop_colname="Proportion_name"){
    
    # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
    Cell_type <- Percentage <- Percentage_label <- NULL
    
    cellprop_colname <- rlang::ensym(cellprop_colname)
    
    cell_proportions$Percentage_label <- 
        round(cell_proportions$Percentage, digits=1)
    cell_proportions <- cell_proportions[cell_proportions$Cell_type != "OTHER",]
    
    cell_percentages_full_plot <-
        ggplot(cell_proportions,
               aes(x = stats::reorder({{cellprop_colname}}, Percentage), 
                   y = Percentage, fill = Cell_type)) +
        geom_bar(stat = 'identity') +
        theme_bw() +
        theme() +
        xlab("Cell Type") + ylab("Proportion of cells") +
        geom_text(aes(label = Percentage_label), 
                  position = position_stack(vjust = 0.5),size = 3) +
        coord_flip()
    
    if(!is.null(cells_to_exclude)){
        cell_proportions_no_tumour <- 
            cell_proportions[!cell_proportions$Cell_type %in% cells_to_exclude,]
        cell_proportions_no_tumour$Proportion <- NULL
        cell_proportions_no_tumour$Percentage <- 
            (cell_proportions_no_tumour$Number_of_celltype/
                 sum(cell_proportions_no_tumour$Number_of_celltype))*100
        cell_proportions_no_tumour$Percentage_label <- 
            round(cell_proportions_no_tumour$Percentage, digits=1)
        
        cell_percentages_no_tumour_plot <-
            ggplot(cell_proportions_no_tumour, 
                   aes(x = stats::reorder({{cellprop_colname}}, Percentage), 
                       y = Percentage, fill=Cell_type)) +
            geom_bar(stat = 'identity') +
            ggtitle(paste("Excluding cells:",cells_to_exclude))+
            theme_bw() +
            theme() +
            xlab("Cell Type") + ylab("Proportion of cells") +
            geom_text(aes(label = Percentage_label), 
                      position = position_stack(vjust = 0.5), size = 2) +
            coord_flip()
        methods::show(cell_percentages_no_tumour_plot)
    }
    methods::show(cell_percentages_full_plot)
    return(cell_percentages_full_plot)
}
