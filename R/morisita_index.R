#' morisita_index (not compatibale with current version; needs update)
#'
#' @description Calculates the morisita index between different formatted single
#'   cell experiment images.
#'
#' @param sce_objects A vector of SingleCellExperiment object names in the form
#'   of the output of format_image_to_sce.
#' @param CI Confidence interval for calculating morisita index.
#' @import dplyr
#' @examples 
#' morisita_index(SPIAT::defined_image)
#' @export

morisita_index <- function(sce_objects, CI = 0.95) {
    
    num_objects <- length(sce_objects)
    
    df_names <- vector()
    
    phenotypes <- vector()

    #FORMAT ALL SCE_OBJECTS
    for(num in 1:num_objects) {
        
        sce_object <- sce_objects[num]
        
        data <- bind_colData_intensity(eval(parse(text = sce_object)))
        
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
    result <- divo::mh(t(count_df), CI=CI)
    
    #set the names of columns and rows
    colnames(result[["Mean"]]) <- sce_objects
    rownames(result[["Mean"]]) <- sce_objects
    colnames(result[["Lower.Quantile"]]) <- sce_objects
    rownames(result[["Lower.Quantile"]]) <- sce_objects    
    colnames(result[["Upper.Quantile"]]) <- sce_objects
    rownames(result[["Upper.Quantile"]]) <- sce_objects
    
    return(result)
}





