#' plot_cell_proportions
#'
#' @description Plots cells proportions as barplots
#' @param cell_proportions Output from calculate_cell_proportions
#' @param tumour_marker Tumour marker to exclude if needed
#' @import ggplot2
#' @export


plot_cell_percentages <- function(cell_proportions, tumour_marker=NULL){
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Cell_type <- Percentage <- Percentage_label <- NULL

  cell_proportions$Percentage_label <- round(cell_proportions$Percentage, digits=1)
  cell_proportions <- cell_proportions[cell_proportions$Cell_type != "OTHER",]

  cell_percentages_full_plot <-
    ggplot(cell_proportions,aes(x = Cell_type,y = Percentage, fill = Cell_type)) +
      geom_bar(stat = 'identity') +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none") +
      xlab("Cell Type") + ylab("Proportion of total cells") +
      geom_text(aes(label = Percentage_label), position=position_dodge(width=0.9), vjust=-0.25, size = 2)

  if(!is.null(tumour_marker)){
    cell_proportions_no_tumour <- cell_proportions[cell_proportions$Cell_type != tumour_marker,]
    cell_proportions_no_tumour$Proportion <- NULL
    cell_proportions_no_tumour$Percentage <- (cell_proportions_no_tumour$Number_of_cells/sum(cell_proportions_no_tumour$Number_of_cells))*100
    cell_proportions_no_tumour$Percentage_label <- round(cell_proportions_no_tumour$Percentage, digits=1)

    cell_percentages_no_tumour_plot <-
      ggplot(cell_proportions_no_tumour, aes(x = Cell_type,y = Percentage,
                                             fill = Cell_type)) +
      geom_bar(stat = 'identity') +
      ggtitle("Excluding tumor cells")+
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none") +
      xlab("Cell Type") + ylab("Proportion of total cells") +
      geom_text(aes(label = Percentage_label), position=position_dodge(width=0.9), vjust=-0.25, size = 2)
    print(cell_percentages_no_tumour_plot)
  }
  print(cell_percentages_full_plot)
}


