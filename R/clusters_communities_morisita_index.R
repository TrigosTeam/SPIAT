#' clusters_communities_morisita_index
#'
#' @description Calculate the morisita index between different clusters or communities
#'
#' @param formatted_data_with_clusters_or_communities A dataframe output from identify_cell_clusters or identify_cell_communities
#' @param column_to_consider Column name to consider as community/clusters
#' @param clusters_or_communities_of_interest A vector containing the names of clusters/communities to calculate morisita index
#' @param CI Confidence interval for calculating morisita index
#' @import dplyr
#' @importFrom divo mh
#' @export

#formatted_data_with_clusters_or_communities <- clusters_data_panimmune_1
#column_to_consider <- "Cluster"
#clusters_or_communities_of_interest <- NULL

#formatted_data_with_clusters_or_communities <- communities_data_panimmune_1
#column_to_consider <- "Community"
#clusters_or_communities_of_interest <- c("Community_1", "Community_2", "Community_7", "Community_4")


clusters_communities_morisita_index <- function(formatted_data_with_clusters_or_communities, column_to_consider, clusters_or_communities_of_interest = NULL, CI = 0.95) {
    
    #SELECT COMMUNITIES TO CALCULATE INDEX IF SPECIFIED
    if (!is.null(clusters_or_communities_of_interest)) {
        formatted_data_with_clusters_or_communities <- formatted_data_with_clusters_or_communities[formatted_data_with_clusters_or_communities$Community %in% clusters_or_communities_of_interest,]
    }
    if(column_to_consider == "Community" || column_to_consider == "Cluster"){
        
        #A 'group' refers to either a community or cluster
        if(column_to_consider == "Community") {
            group_names <- unique(formatted_data_with_clusters_or_communities$Community)
        } else {
            group_names <- unique(formatted_data_with_clusters_or_communities$Cluster)
        }
        num_groups <- length(group_names)
        phenotypes <- unique(formatted_data_with_clusters_or_communities$Phenotype)
        
        #INITIALIZE A DATAFRAME TO STORE THE COUNTS
        count_df <- data.frame(matrix(nrow = num_groups, ncol = length(phenotypes)))
        colnames(count_df) <- phenotypes
        
        #Rename the cluster/community column to 'group'
        colnames(formatted_data_with_clusters_or_communities)[12] <- "Group"
        
        #READ IN FREQ, row order in count_df is the order in group_names
        rownum = 1
        for (group_name in group_names) {
            
            group_data <- formatted_data_with_clusters_or_communities[formatted_data_with_clusters_or_communities$Group == group_name, ]
            
            for (phenotype in phenotypes) {
                #print(phenotype)
                
                if (!(phenotype %in% group_data$Phenotype)) {
                    count <- 0
                } else {
                    count <- nrow(group_data[group_data$Phenotype == phenotype,])
                }
                #print(count)
                count_df[rownum,phenotype] <- count
            }
            rownum = rownum + 1
        }
        
        #START A NEW DATAFRAME FOR STORING THE MORISITA INDEX
        result_df <- data.frame(matrix(data = NA, nrow=num_groups, ncol=num_groups))
        colnames(result_df) <- group_names
        rownames(result_df) <- group_names
        
        count_df_t <- t(count_df)
        
        #TRY CATCH COMMUNITIES WITH LITTLE CELLS THAT LEADS TO ERROR
        calc_morisita <- function(matrix_of_two_groups, CI, group1_name, group2_name) {
            
            out <- tryCatch(
                {
                    mh(matrix_of_two_groups, CI=CI)
                }, 
                error=function(cond) {
                    message(paste("Cannot calculate morisita index between ", group1_name, " and ", group2_name, sep=""))
                    return(NULL)
                }, 
                warning=function(cond) {
                    return(NULL)
                }, 
                finally={
                }
            )
            return(out)
        }
        
        #PAIRWISE CALCULATION OF MORISITA INDEX AND UPDATE RESULT
        for (group1 in 1:num_groups) {
            group1_name <- group_names[group1]
            
            for (group2 in 1:num_groups) {
                if (group2 < group1) {
                    next
                } else if (group2 == group1) {
                    result_df[group2, group1] <- 1
                    next
                }
                group2_name <- group_names[group2]
                
                temp <- count_df_t[,c(group1,group2)]
                
                out <- calc_morisita(temp, CI=CI, group1_name, group2_name)
                
                if (!is.null(out)) {
                    mean_index <- out$Mean[2,1]
                    result_df[group2, group1] <- mean_index
                }
                
            }
        }
        
        #result <- mh(t(count_df), CI=CI)
        
        #set the names of columns and rows
        #colnames(result[["Mean"]]) <- group_names
        #rownames(result[["Mean"]]) <- group_names
        #colnames(result[["Lower.Quantile"]]) <- group_names
        #rownames(result[["Lower.Quantile"]]) <- group_names    
        #colnames(result[["Upper.Quantile"]]) <- group_names
        #rownames(result[["Upper.Quantile"]]) <- group_names
        
        return(result_df)
        
    }else{
        stop("Only Community and Cluster are accepted as valid column names")
    }
       
}