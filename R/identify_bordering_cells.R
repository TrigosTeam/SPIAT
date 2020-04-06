#' identify_bordering_cells
#'
#' @description Identifies the cells bordering a group of cells of a particular phenotype
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param marker Cells positive for this marker will be used as reference
#' @param rm_noise_radius Number specifying the radius used for noise. Larger numbers, more leniency
#' @param radius Number specifying the search radius. Larger numbers, more cells to be considered
#' @param lower_bound Number specifying the minumum proportion of non-marker cells in the radius of the marker population
#' @param upper_bound Number specifying the maximum proportion of non-marker cells in the radius of the marker population
#' @import SingleCellExperiment
#' @export

#sce_object <- sce_ovarian_panimmune1
#tumour_marker <- "WT1"
#rm_noise_radius <- 50 #larger the more lenient it is
#radius <- 200 #larger the more cells are considered
#lower_bound <- 0.4
#upper_bound <- 0.6

identify_bordering_cells <- function(sce_object, tumour_marker, rm_noise_radius, radius, lower_bound, upper_bound) {

  formatted_data <- data.frame(colData(sce_object))

  formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

  expression_matrix <- assay(sce_object)

  markers <- rownames(expression_matrix)
  cell_ids <- colnames(expression_matrix)

  rownames(expression_matrix) <- NULL
  colnames(expression_matrix) <- NULL
  expression_matrix_t <- t(expression_matrix)
  expression_df <- data.frame(expression_matrix_t)
  colnames(expression_df) <- markers

  formatted_data <- cbind(formatted_data, expression_df)
  formatted_data <- formatted_data[complete.cases(formatted_data),]

  #####################

  #get the tumour cells and other cells
  tumour_cells <- formatted_data[formatted_data$Phenotype == tumour_marker, ]
  other_cells <- formatted_data[!grepl(tumour_marker, formatted_data$Phenotype), ]

  #remove tumour cells with no tumour cells within the radius
  tumour_cell_cords <- tumour_cells[,c("Cell.X.Position", "Cell.Y.Position")]
  tumour_nn <- frNN(tumour_cell_cords, eps=rm_noise_radius, sort=FALSE)  #row index matching the row num
  #get the row index of those cells with neighbours
  idx <- unique(unlist(tumour_nn$id))
  tumour_cells_w_neighbours <- tumour_cells[idx, ]

  #remove stroma cells with no stroma cells within a specific radius
  other_cell_cords <- other_cells[,c("Cell.X.Position", "Cell.Y.Position")]
  other_nn <- frNN(other_cell_cords, eps=rm_noise_radius, sort=FALSE)
  idx <- unique(unlist(other_nn$id))
  other_cells_w_neighbours <- other_cells[idx, ]

  #create df of all cells
  all_cells <- rbind(tumour_cells_w_neighbours, other_cells_w_neighbours)
  all_cell_cords <- all_cells[,c("Cell.X.Position", "Cell.Y.Position")]

  #use tumour cells as reference, search for neighbours and get tumour cells with 40~60% stroma cells around
  all_nn <- frNN(x=all_cell_cords, eps=radius, query=tumour_cell_cords, sort=FALSE)
  tumour_row_nums <- rownames(tumour_cells)
  border_cells <- data.frame(matrix(ncol = ncol(tumour_cells), nrow = 0))
  colnames(border_cells) <- colnames(tumour_cells)
  for (row_num in tumour_row_nums) {
    nn_index <- all_nn$id[row_num]
    unlisted_index <- unique(unlist(nn_index))
    nn_num = length(unlisted_index)
    if (nn_num <= 1){
      next
    }

    nn <- all_cells[unlisted_index, ]

    non_tumour_count <- nrow(nn[!grepl(tumour_marker, nn$Phenotype), ])

    prop_non_tumour <- non_tumour_count/nn_num

    if (prop_non_tumour >= lower_bound && prop_non_tumour <= upper_bound) {
      border_cells <- rbind(border_cells, tumour_cells[row_num, ])
    }
  }

  r <- ggplot(border_cells, aes(x = Cell.X.Position, y = Cell.Y.Position)) +
    geom_point(size = 0.1) +
    guides(alpha = F) + scale_colour_viridis_c(direction = -1) +
    labs(colour = paste("log10","(", as.character(tumour_marker)," Intensity", ")", sep="")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), legend.key.height = unit(2.5, "cm"))
  print(r)
}
