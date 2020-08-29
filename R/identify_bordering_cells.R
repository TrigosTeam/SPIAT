#' identify_bordering_cells
#'
#' @description Identifies the cells bordering a group of cells of a particular phenotype
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param reference_marker Cells positive for this marker will be used as reference
#' @param rm_noise_radius Number specifying the radius used for noise. Larger numbers, more leniency
#' @param radius Number specifying the search radius. Larger numbers, more cells to be considered
#' @param lower_bound Number specifying the minumum proportion of non-marker cells in the radius of the reference marker population
#' @param upper_bound Number specifying the maximum proportion of non-marker cells in the radius of the reference marker population
#' @importFrom SummarizedExperiment assay colData
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom dbscan frNN
#' @import ggplot2
#' @return A plot is returned
#' @examples
#' identify_bordering_cells(SPIAT::formatted_image, reference_marker = "AMACR",
#'                          rm_noise_radius = 50, radius = 100, lower_bound = 0.05,
#'                          upper_bound=0.7)
#' @export

#sce_object <- sce_ovarian_panimmune1
#reference_marker <- "WT1"
#rm_noise_radius <- 50 #larger the more lenient it is
#radius <- 200 #larger the more cells are considered
#lower_bound <- 0.4
#upper_bound <- 0.6

# imported all functions from ggplot2 since many functions are interdependent

identify_bordering_cells <- function(sce_object, reference_marker, rm_noise_radius, radius, lower_bound, upper_bound) {
  
  # setting these column names to NULL as otherwise get "no visible binding for global variable" in R check
  Cell.X.Position <- Cell.Y.Position <- NULL

  formatted_data <- data.frame(colData(sce_object))

  formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

  intensity_matrix <- assay(sce_object)

  markers <- rownames(intensity_matrix)
  cell_ids <- colnames(intensity_matrix)

  rownames(intensity_matrix) <- NULL
  colnames(intensity_matrix) <- NULL
  intensity_matrix_t <- t(intensity_matrix)
  intensity_df <- data.frame(intensity_matrix_t)
  colnames(intensity_df) <- markers

  formatted_data <- cbind(formatted_data, intensity_df)
  formatted_data <- formatted_data[complete.cases(formatted_data),]

  #####################

  #get the reference cells and other cells
  reference_cells <- formatted_data[formatted_data$Phenotype == reference_marker, ]
  other_cells <- formatted_data[!grepl(reference_marker, formatted_data$Phenotype), ]
  
  #CHECK
  if (nrow(reference_cells) == 0){
    stop("There are no reference_cells in the dataset for the specified marker")
  }
  if (nrow(other_cells) == 0){
    stop("There are no stroma cells in the dataset")
  }
  
  #remove reference cells with no reference cells within the radius
  reference_cell_cords <- reference_cells[,c("Cell.X.Position", "Cell.Y.Position")]
  reference_nn <- frNN(reference_cell_cords, eps=rm_noise_radius, sort=FALSE)  #row index matching the row num
  #get the row index of those cells with neighbours
  idx <- unique(unlist(reference_nn$id))
  reference_cells_w_neighbours <- reference_cells[idx, ]
  #CHECK
  if (nrow(reference_cells_w_neighbours) == 0) {
    stop("No reference cells left after noise removal, please try a more lenient radius")
  }

  #remove stroma cells with no stroma cells within a specific radius
  other_cell_cords <- other_cells[,c("Cell.X.Position", "Cell.Y.Position")]
  other_nn <- frNN(other_cell_cords, eps=rm_noise_radius, sort=FALSE)
  idx <- unique(unlist(other_nn$id))
  other_cells_w_neighbours <- other_cells[idx, ]
  #CHECK
  if (nrow(other_cells_w_neighbours) == 0) {
    stop("No stroma cells left after noise removal, please try a more lenient radius")
  }

  #create df of all cells
  all_cells <- rbind(reference_cells_w_neighbours, other_cells_w_neighbours)
  all_cell_cords <- all_cells[,c("Cell.X.Position", "Cell.Y.Position")]

  #use reference_cells, search for neighbours and get reference cells with 40~60% stroma cells around
  all_nn <- frNN(x=all_cell_cords, eps=radius, query=reference_cell_cords, sort=FALSE)
  reference_row_nums <- rownames(reference_cells)
  border_cells <- data.frame(matrix(ncol = ncol(reference_cells), nrow = 0))
  colnames(border_cells) <- colnames(reference_cells)
  for (row_num in reference_row_nums) {
    nn_index <- all_nn$id[row_num]
    unlisted_index <- unique(unlist(nn_index))
    nn_num = length(unlisted_index)
    if (nn_num <= 1){
      next
    }

    nn <- all_cells[unlisted_index, ]

    non_reference_count <- nrow(nn[!grepl(reference_marker, nn$Phenotype), ])

    prop_non_reference <- non_reference_count/nn_num

    if (prop_non_reference >= lower_bound && prop_non_reference <= upper_bound) {
      border_cells <- rbind(border_cells, reference_cells[row_num, ])
    }
  }
  
  #CHECK
  if (nrow(border_cells) == 0){
    stop("There are no border cells found for the specified radius")
  }

  r <- ggplot(border_cells, aes(x = Cell.X.Position, y = Cell.Y.Position)) +
    geom_point(size = 0.1) +
    guides(alpha = FALSE) + scale_colour_viridis_c(direction = -1) +
    labs(colour = paste("log10","(", as.character(reference_marker)," Intensity", ")", sep="")) +
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
