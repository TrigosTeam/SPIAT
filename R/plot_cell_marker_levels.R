#' plot_cell_marker_levels
#'
#' @description Produces a scatter plot of the expression of every marker in each cell.
#' Cells that were not phenotyped as being positive for the particular marker are excluded.
#'
#' @param sce_object Singlecellexperiment object in the form of the output of format_image_to_sce
#' @param print TRUE for plots to be printed, FALSE otherwise
#' @param filename Path and name of output pdf file if to be printing to a file, if required
#' @param return_data TRUE if the function should return the formatted data for plotting
#' @import SingleCellExperiment
#' @import qpdf
#' @import ggplot2
#' @import grDevices
#' @export

plot_cell_marker_levels <- function(sce_object, print = TRUE, filename=NULL, return_data=TRUE) {

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

    for (marker in markers) {

        #every cell is stained by DAPI, so no need to remove intensity
        if (marker == "DAPI"){
            next
        }

        #selecting cells that do not contain the marker
        rows <- formatted_data[formatted_data$Phenotype != marker, ] #for one entry that is not marker
        rows <- rows[!grepl(marker, rows$Phenotype), ] #for multiple entries that does not contain marker

        #for those cell without the marker, set marker intensity to 0
        #and merge the formatted_data
        rows[, marker] <- 0
        formatted_data[match(rows$Cell.ID,formatted_data$Cell.ID),]<-rows

    }



    if(print){
      if(!is.null(filename)){
        pdf(filename)
      }

      for (marker in markers){

        #selecting the cells that have intensity for a specific marker
        column <- which(colnames(formatted_data) == marker)
        rows_non_zero <- which(formatted_data[,column] != 0)
        intensity_by_marker <- formatted_data[rows_non_zero,]

        if (nrow(intensity_by_marker) == 0) {
          print(paste("There are no true expression for: ", marker, sep=""))
        }

        #log the intensity to improve contrast
        intensity_by_marker[,marker] <- log10(intensity_by_marker[,marker])
        #print(intensity_by_marker)

        p <- ggplot(intensity_by_marker, aes(x = Cell.X.Position, y = Cell.Y.Position, colour = eval(parse(text = marker)))) +
          geom_point(aes(colour=eval(parse(text = marker))),size = 0.1) +
          ggtitle(marker)+
          guides(alpha = F) + scale_colour_viridis_c(direction = -1) +
          labs(colour = paste("log10","(", as.character(marker)," Intensity", ")", sep="")) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "white"),
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(), legend.key.height = unit(2.5, "cm"))

        print(p)
      }
    }

    if(!is.null(filename)){
      dev.off(filename)
    }

    if(isTRUE(return_data)){
      return(formatted_data)
    }
}

