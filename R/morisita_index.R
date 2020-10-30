#' morisita_index
#'
#' @description Calculate the morisita index between different formatted single cell experiment images
#'
#' @param sce_objects A vector of SingleCellExperiment object names in the form of the output of format_image_to_sce
#' @param CI Confidence interval for calculating morisita index
#' @import dplyr
#' @import SingleCellExperiment
#' @importFrom tibble rownames_to_column
#' @importFrom divo mh
#' @export

morisita_index <- function(sce_objects, CI = 0.95) {
    
    num_objects <- length(sce_objects)
    
    df_names <- vector()
    
    phenotypes <- vector()
    
    format_sce_to_df <- function(sce_object) {
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
        
        return(formatted_data)
    }
    
    #FORMAT ALL SCE_OBJECTS
    for(num in 1:num_objects) {
        
        sce_object <- sce_objects[num]
        
        data <- format_sce_to_df(eval(parse(text = sce_object)))
        
        phenotypes <- unique(c(phenotypes, unique(data$Phenotype)))
        
        assign(paste("formatted_data", "_", as.character(num), sep=""), data)
        
        df_names <- c(df_names, paste("formatted_data", "_", as.character(num), sep=""))
    }
    
    #INITIALIZE A DATAFRAME TO STORE THE COUNTS
    count_df <- data.frame(matrix(nrow = num_objects, ncol = length(phenotypes)))
    colnames(count_df) <- phenotypes
    
    #POPULATE FREQUENCY
    for (num in 1:num_objects) {
        
        df <- data.frame(get(df_names[num]))
        
        for (phenotype in phenotypes) {
            #print(phenotype)
            
            if (!(phenotype %in% df$Phenotype)) {
                count <- 0
            } else {
                count <- nrow(df[df$Phenotype == phenotype,])
            }
            #print(count)
            count_df[num,phenotype] <- count
        }
    }
    
    #CALCULATE MORISITA INDEX
    result <- mh(t(count_df), CI=CI)
    
    #set the names of columns and rows
    colnames(result[["Mean"]]) <- sce_objects
    rownames(result[["Mean"]]) <- sce_objects
    colnames(result[["Lower.Quantile"]]) <- sce_objects
    rownames(result[["Lower.Quantile"]]) <- sce_objects    
    colnames(result[["Upper.Quantile"]]) <- sce_objects
    rownames(result[["Upper.Quantile"]]) <- sce_objects
    
    return(result)
}





