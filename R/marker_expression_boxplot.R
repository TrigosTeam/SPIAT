#' marker_expression_boxplot
#'
#' @description Produces boxplots of marker levels for cells phenotyped as being positive
#' for the marker, and those that where phenotyped as being negative.
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param marker Marker being queried
#' @importFrom SummarizedExperiment colData assay
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom stats complete.cases
#' @import ggplot2
#' @export

marker_expression_boxplot <- function(sce_object, marker){

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

    #selecting cells that do express the marker
    expression_true <- formatted_data[grepl(marker, formatted_data$Phenotype), ]

    #selecting cells that do not contain the marker
    expression_false <- formatted_data[!grepl(marker, formatted_data$Phenotype), ] #for multiple entries that does not contain marker

    #select the specific marker and add a boolean of expression
    if (nrow(expression_true) != 0) {
        expression_true$expression <- "P"
    } else{
        stop(paste("There are no cells positive for ", marker, sep=""))
    }

    if (nrow(expression_false) != 0){
        expression_false$expression <- "N"
    } else{
        stop(paste("There are no cells negative for ", marker, sep=""))
    }

    #bind the 2 dataframes together
    expression_by_marker <- rbind(expression_true, expression_false)
    expression_by_marker$expression <- factor(expression_by_marker$expression, levels=c("P", "N"))
    
    #plot boxplot
    title <- paste("Level of ", marker, sep="")

    #code from ggplot2 to show the number in each group
    give.n <- function(x){
        return(c(y = max(x)+1, label = length(x)))
    }

    p <- ggplot(expression_by_marker, aes(x = expression, y = expression_by_marker[, marker]))
    p <- p + geom_boxplot()
    p <- p + labs(title = title, x = "Marker status (Positive/Negative cell)", y = "Marker level")
    p <- p + stat_summary(fun.data = give.n, geom = "text")
    p <- p + theme_bw()

    print(p)
}
