#' predict_phenotypes
#'
#' @description Produces a density plot showing actual and predicted cutoff of a
#' positive reading for marker expression. It also prints to the console of the
#' number of true positives (TP), true negatives (TN), false positives (FP) and
#' false negatives (FN) under the prediction. It returns a dataframe containing
#' the predicted expression status for a particular marker,

#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param thresholds (Optional) Vector of numbers specifying the cutoff of a positive reading.
#' The order must match the marker order, and it should be NA for DAPI.
#' @param tumour_marker String containing the tumour_marker used for the image.
#' @param baseline_markers Markers not found on tumour cells to refine the threshold
#' @import SingleCellExperiment
#' @import pracma
#' @export

predict_phenotypes <- function(sce_object, plot_actual_cutoff = FALSE, plot_predicted_cutoff = FALSE, thresholds = NULL, tumour_marker,
                               baseline_markers) {

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

    #add actual expression boolean value to formatted_data

    markers_no_tumour <- markers[markers != tumour_marker]

    ##Add actual marker levels
    for (marker in markers){
        if (marker == "DAPI"){
            actual_phenotype <- data.frame(rep(1, nrow(formatted_data)))
            colnames(actual_phenotype) <- paste(marker, "_actual_phenotype", sep="")
            formatted_data <- cbind(formatted_data, actual_phenotype)
        } else{
            cell_IDs <- data.frame(formatted_data$Cell.ID)
            colnames(cell_IDs) <- "Cell.ID"
            actual_phenotype <- data.frame(rep(0, nrow(formatted_data)))
            colnames(actual_phenotype) <- "Phenotype_status"
            actual_phenotype <- cbind(actual_phenotype, cell_IDs)

            rows <- formatted_data[grepl(marker, formatted_data$Phenotype), ]
            actual_phenotype[actual_phenotype$Cell.ID %in% rows$Cell.ID,]$Phenotype_status <- 1

            actual_phenotype <- data.frame(actual_phenotype[,1])
            colnames(actual_phenotype) <- paste(marker, "_actual_phenotype", sep="")
            formatted_data <- cbind(formatted_data, actual_phenotype)
        }
    }

    for (marker in markers_no_tumour) {
        if (marker == "DAPI") {
            next
        }
        #extract the marker expression column
        marker_specific_level <- formatted_data[,marker]

        #calculate the predictions
        if (!is.null(thresholds) && !is.na(thresholds[match(marker,markers)])) {
            #there is a threshold value specified for the marker, use the threshold
            marker_threshold <- thresholds[match(marker,markers)]
            print(paste("(", marker, " has threshold specified: ", as.character(marker_threshold), ")", sep=""))
            selected_valley_xcord <- NULL

            #get the threshold predictions
            predictions_by_threshold <- data.frame(threshold(marker_specific_level, level = marker_threshold))

        } else {
            #calculate the valleys
            expression_density <- density(marker_specific_level)
            valleys <- findpeaks(-(expression_density)$y)
            valley_ycords <- valleys[,1] * -1
            index <- match(valley_ycords, expression_density$y)
            valley_xcords <- expression_density$x[index]

            #create a df for the valley coordinates
            valley_df <- data.frame(cbind(valley_xcords, valley_ycords))

            #select the first valley that's greater than the maximum density and below 25% density
            ycord_max_density <- max(expression_density$y)
            xcord_max_density_index <- match(ycord_max_density, expression_density$y)
            xcord_max_density <- expression_density$x[xcord_max_density_index]

            density_threshold_for_valley <- 0.25 * ycord_max_density

            valley_df <- valley_df[valley_df$valley_xcords >= xcord_max_density, ]
            valley_df <- valley_df[valley_df$valley_ycords <= density_threshold_for_valley, ]

            selected_valley_xcord <- valley_df$valley_xcords[1]
            selected_valley_ycord <- valley_df$valley_ycords[1]
        }
            #using the selected valley as the threshold
            predictions_by_threshold <- data.frame(threshold(marker_specific_level, level = selected_valley_xcord))
            colnames(predictions_by_threshold) <- paste(marker, "_predicted_phenotype", sep="")
            formatted_data <- cbind(formatted_data, predictions_by_threshold)
      }


    ###Prediction for tumour marker
    #Select cells that are positive for a marker that tumor cells don't have

    #baseline_markers <- c("FOXP3", "CD3", "CD4", "CD8")


    baseline_cells <- vector()
    for(marker in baseline_markers){
      temp <- formatted_data[, colnames(formatted_data) %in% c("Cell.ID", paste0(marker, "_predicted_phenotype"))]
      temp <- temp[temp[,2] == 1,]
      baseline_cells <- c(baseline_cells, temp$Cell.ID)
    }
    baseline_cells <- unique(baseline_cells)

    #Tumor marker levels in these cells

    formatted_data_baseline <- formatted_data[formatted_data$Cell.ID %in% baseline_cells,tumour_marker]
    cutoff_for_tumour <- quantile(formatted_data_baseline, 0.95)

    #extract the marker expression column
    tumour_specific_level <- formatted_data[,tumour_marker]

    #calculate the predictions
    if (!is.null(thresholds)) {
      #there is a threshold value specified for the marker, use the threshold
      marker_threshold <- thresholds[match(tumour_,markers)]
      print(paste("(", tumour_marker, " has threshold specified: ", as.character(marker_threshold), ")", sep=""))
      selected_valley_xcord <- NULL

      #get the threshold predictions
      predictions_by_threshold <- data.frame(threshold(tumour_specific_level, level = marker_threshold))

    } else {
      #calculate the valleys
      expression_density <- density(tumour_specific_level)
      valleys <- findpeaks(-(expression_density)$y)
      valley_ycords <- valleys[,1] * -1
      index <- match(valley_ycords, expression_density$y)
      valley_xcords <- expression_density$x[index]

      #create a df for the valley coordinates
      valley_df <- data.frame(cbind(valley_xcords, valley_ycords))
      selected_valley_xcord <- valley_df$valley_xcords[1]
      selected_valley_ycord <- valley_df$valley_ycords[1]

      #using the selected valley as the threshold if it is lower than the
      #level of expression of the tumour marker in non-tumour cells

      final_threshold <- ifelse(selected_valley_xcord < cutoff_for_tumour,
                                selected_valley_xcord, cutoff_for_tumour)
      predictions_by_threshold <- data.frame(threshold(tumour_specific_level, level = final_threshold))
      colnames(predictions_by_threshold) <- paste(tumour_marker, "_predicted_phenotype", sep="")
      formatted_data <- cbind(formatted_data, predictions_by_threshold)
    }
    predicted_data <- formatted_data
    
    markers <- markers[markers != "DAPI"]
    
    for(marker in markers){
      #names of columns
      marker_status_name <- paste(marker, "_status", sep="")
      marker_actual_exp_colname <- paste(marker,"_actual_phenotype", sep="")
      marker_pred_exp_colname <- paste(marker,"_predicted_phenotype", sep="")
      
      #create an accuracy_df with the same number of rows and 1 column
      accuracy_df <- data.frame(rep(NA, nrow(predicted_data)))
      colnames(accuracy_df) <- "status"
      
      #grab both the actual and predicted expression for the specific marker, change colnames and bind to accuracy_df
      marker_exp_actual_pred <- predicted_data[,c(marker_actual_exp_colname, marker_pred_exp_colname)]
      colnames(marker_exp_actual_pred) <- c("actual", "pred")
      accuracy_df <- cbind(accuracy_df, marker_exp_actual_pred)
      
      #set the TP, TF, FP, FN if they exist
      if (nrow(accuracy_df[accuracy_df$actual == 1 & accuracy_df$pred == 1 , ])) {
        accuracy_df[accuracy_df$actual == 1 & accuracy_df$pred == 1 , ]$status <- "TP"
      }
      if (nrow(accuracy_df[accuracy_df$actual == 0 & accuracy_df$pred == 0 , ])) {
        accuracy_df[accuracy_df$actual == 0 & accuracy_df$pred == 0 , ]$status <- "TN"
      }
      if (nrow(accuracy_df[accuracy_df$actual == 0 & accuracy_df$pred == 1 , ])) {
        accuracy_df[accuracy_df$actual == 0 & accuracy_df$pred == 1 , ]$status <- "FP"
      }
      if (nrow(accuracy_df[accuracy_df$actual == 1 & accuracy_df$pred == 0 , ])) {
        accuracy_df[accuracy_df$actual == 1 & accuracy_df$pred == 0 , ]$status <- "FN"
      }
      
      #bind the specific marker_status to expression level
      accuracy_df <- data.frame(accuracy_df[,"status"])
      colnames(accuracy_df) <- "status"
      marker_specific_level <- predicted_data[,marker]
      marker_specific_level <- data.frame(marker_specific_level)
      colnames(marker_specific_level) <- "Marker_level"
      level_and_accuracy <- cbind(marker_specific_level, accuracy_df)
      
      #print the number of TP, TN, FP, FN
      TP_count <- nrow(level_and_accuracy[level_and_accuracy$status == "TP", ])
      TN_count <- nrow(level_and_accuracy[level_and_accuracy$status == "TN", ])
      FP_count <- nrow(level_and_accuracy[level_and_accuracy$status == "FP", ])
      FN_count <- nrow(level_and_accuracy[level_and_accuracy$status == "FN", ])
      print(paste("For ", marker, ":", sep=""))
      print(paste("TP:", TP_count, " TN:", TN_count, " FP:", FP_count, " FN:", FN_count, sep=""))
      
      #determine the threshold
      real_pos <- level_and_accuracy[level_and_accuracy$status %in% c("TP", "FN"), ]
      threshold_val_real_pos <- min(real_pos$Marker_level)
      
      real_neg <- level_and_accuracy[level_and_accuracy$status %in% c("FP", "TN"), ]
      threshold_val_real_neg <- max(real_neg$Marker_level)
      
      p <- ggplot(level_and_accuracy, aes(x=Marker_level)) + geom_density()
      title <- paste("Density distribution of", marker, sep=" ")
      p <- p + labs(title = title, x = "Level of expression", y = "Density")
      
      if (plot_actual_cutoff == TRUE) {
        p <- p + geom_vline(aes(xintercept=threshold_val_real_pos, color="min_pos"), alpha = 0.3)
        p <- p + geom_vline(aes(xintercept=threshold_val_real_neg, color="max_neg"), alpha = 0.3)
      }
      if (plot_predicted_cutoff == TRUE) {
        if (!is.null(selected_valley_xcord)) {
          p <- p + geom_vline(aes(xintercept = selected_valley_xcord, color = "cutoff"), alpha = 0.3)
        } else {
          p <- p + geom_vline(aes(xintercept = marker_threshold, color = "cutoff"), alpha = 0.3)
        }
      }
      p <- p + scale_color_manual(name = "lines", values = c(min_pos = "red", max_neg = "blue", cutoff = "green"))
      
      print(p)
    }
    return (predicted_data)
}
