#' marker_permutation
#'
#' @description Creates random combinations of phenotypes by shuffling markers and
#' calculates the enrichment and depletion p values
#' @param sce_object SingleCellExperiment object in the form of output from format_image_to_sce
#' @param num_iter Integer specifying the number of iterations for bootstrapping
#' @import dplyr
#' @importFrom SummarizedExperiment colData assay
#' @importFrom tibble rownames_to_column
#' @importFrom stats complete.cases
#' @importFrom utils combn
#' @return A plot is returned
#' @examples 
#' sig <- marker_permutation(SPIAT::formatted_image, num_iter = 100)
#' @export

marker_permutation <- function(sce_object, num_iter) {
  
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
  
  markers <- markers[markers != "DAPI"]
  
  #generate all combinations of markers into a vector
  marker_combinations <- vector()
  for(i in seq_along(markers)) {
    comb_matrix <- combn(markers, i)
    for(j in seq_len(ncol(comb_matrix))) {
      comb <- paste0(comb_matrix[,j], collapse = '', sep=',')
      comb <- gsub(",$", "", comb)
      marker_combinations <- c(marker_combinations, comb)
    }
  }
  
  #create the results df to store the output of every bootstrap iteration
  results <- data.frame(matrix(0, nrow = length(marker_combinations), ncol = num_iter))
  rownames(results) <- marker_combinations
  colnames(results) <- seq_len(num_iter)
  
  #count the markers and put counts into a df
  marker_count <- vector()
  for (marker in markers) {
    count <- nrow(formatted_data[grepl(marker, formatted_data$Phenotype), ])
    marker_count <- c(marker_count, count)
  }
  count_df <- data.frame(matrix(marker_count, ncol=length(marker_count), nrow = 1))
  colnames(count_df) <- markers
  
  total_cells <- nrow(formatted_data)
  
  #ITERATION LOOP here...##############
  for (iter_num in seq_len(num_iter)) {
    #bootstrap_df is used in the randomization of markers to generate random phenotypes
    bootstrap_df <- data.frame(matrix(0, nrow = nrow(formatted_data), ncol = length(markers)))
    rownames(bootstrap_df) <- formatted_data$Cell.ID
    colnames(bootstrap_df) <- markers
    
    #randomly assign marker intensities
    for (marker in markers) {
      marker_count <- count_df[,marker]
      rand <- sample(total_cells, size = marker_count)
      bootstrap_df[rand,marker] <- 1
    }
    
    #start a new column for phenotype
    bootstrap_df$Phenotype <- ""
    for (i in seq_along(markers)) {
      marker <- paste(markers[i], ",", sep="")
      #select the marker column that was assigned to be true (1) for the marker and add the marker name as phenotype
      if(sum(bootstrap_df[, i] == 1) > 0){
        bootstrap_df[bootstrap_df[, i] == 1, ]$Phenotype <- paste(bootstrap_df[bootstrap_df[, i] == 1, ]$Phenotype, marker, sep="")
      }
    }
    
    #get rid of comma at the end
    bootstrap_df$Phenotype <- gsub(",$", "", bootstrap_df$Phenotype)
    #add "OTHER" as phenotype for those without a phenotype, since they're all DAPI positive
    bootstrap_df[bootstrap_df$Phenotype == "", ]$Phenotype <- "OTHER"
    
    #get all unique phenotypes generated
    phenotypes_generated <- unique(bootstrap_df$Phenotype)
    #count the phenotypes and read it into result df
    for (phenotype in phenotypes_generated) {
      count <- nrow(bootstrap_df[bootstrap_df$Phenotype == phenotype, ])
      results[phenotype, iter_num] <- count
    }
    
  }
  
  #start a summary dataframe
  summary_df <- data.frame(matrix(nrow = length(marker_combinations), ncol=5))
  rownames(summary_df) <- marker_combinations
  colnames(summary_df) <- c("Observed_cell_number", "Percentage_of_iterations_where_present", "Average_bootstrap_cell_number",
                            "Enrichment.p", "Depletion.p")
  
  
  #calculate the percentage, enrichment, depletion, observed_cell_number, average_bootstrap_cell_number
  for (combination in marker_combinations) {
    
    #observed_cell_number
    num_observed <- nrow(formatted_data[formatted_data$Phenotype == combination, ])
    summary_df[combination, "Observed_cell_number"] <- num_observed
    
    #percentage
    combination_scores <- results[combination, ]
    num_non_zero <- length(combination_scores[combination_scores != 0])
    percentage_presence <- num_non_zero/num_iter * 100
    summary_df[combination, "Percentage_of_iterations_where_present"] <- percentage_presence
    
    #average_bootstrap_cell_number
    average_bootstrap_cell_number <- mean(unlist(combination_scores))
    summary_df[combination, "Average_bootstrap_cell_number"] <- average_bootstrap_cell_number
    
    #enrichment
    num_greater <- sum(num_observed > combination_scores)
    if (length(num_greater) == 0) {
      enrichment_score <- 1
    } else {
      enrichment_score <- 1-(num_greater/num_iter)
      if(enrichment_score == 0){
        enrichment_score <- 1/num_iter
      }
    }
    summary_df[combination, "Enrichment.p"] <- enrichment_score
    
    #depletion
    num_lesser <- sum(num_observed < combination_scores)
    if (length(num_lesser) == 0) {
      depletion_score <- 1
    } else {
      depletion_score <- 1-(num_lesser/num_iter)
      if(depletion_score == 0){
        depletion_score <- 1/num_iter
      }
    }
    summary_df[combination, "Depletion.p"] <- depletion_score
  }
  
  return(summary_df)
  
}
