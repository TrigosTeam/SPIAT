#' calculate_summary_distances_between_phenotypes
#'
#' @description Returns the mean, median and stardard deviation of the distances between phenotypes
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param all_marker_combinations If TRUE, the distances between all possible combinations of markers
#' will be calculated
#' @param combinations_of_interest Vector of marker combinations to consider if all_marker_combinations
#' is FALSE
#' @importFrom apcluster negDistMat
#' @importFrom RANN nn2
#' @importFrom gtools permutations
#' @import dplyr
#' @importFrom stats median sd
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble rownames_to_column
#' @return A data.frame is returned
#' @examples
#' summary_distances <- calculate_summary_distances_between_phenotypes(SPIAT::formatted_image)
#' @export

calculate_summary_distances_between_phenotypes <- function(sce_object, all_marker_combinations = TRUE,
                                                   combinations_of_interest = NULL) {

    formatted_data <- data.frame(colData(sce_object))
    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

    formatted_data <- formatted_data[,c("Cell.ID","Phenotype", "Cell.X.Position", "Cell.Y.Position")]
    formatted_data <- formatted_data[formatted_data$Phenotype != "",]

    #Add a new column Cell_type which duplicates the phenotype
    formatted_data$Cell_type <- as.character(formatted_data$Phenotype)

    #Get the list of cells under each cell type
    cell_types = list()
    if (all_marker_combinations) {
        #assign cells to each cell type based on the markers expressed
        for (eachType in unique(formatted_data$Cell_type)) {
            cell_types[[eachType]] = as.character(formatted_data$Cell.ID[formatted_data$Phenotype == eachType])
        }
        print("All markers are used in pair-wise distance calculation: ")
        print(unique(formatted_data$Cell_type))
    }else{
        #assign cells to each specified cell type
        for (eachType in combinations_of_interest) {
            cell_types[[eachType]] = as.character(formatted_data$Cell.ID[formatted_data$Phenotype == eachType])
        }
        #keep those cells that are selected
        formatted_data <- formatted_data[formatted_data$Phenotype %in% combinations_of_interest, ]
        #CHECK
        if (nrow(formatted_data) == 0){
            stop("No cells belong to the specified marker combinations of interest")
        }
        print("Markers had been selected in pair-wise distance calculation: ")
        print(unique(formatted_data$Cell_type))
    }

    #permutation for the different cell type combinations
    permu = permutations(length(unique(formatted_data$Cell_type)), 2, repeats.allowed = TRUE)
    unique_cells <- unique(formatted_data$Cell_type) #unique cell types
    result = vector()

    for (i in seq_len(nrow(permu))) {
        eachPermu = permu[i, ]
        name1 = unique_cells[eachPermu[1]]
        name2 = unique_cells[eachPermu[2]]

        cell_type1 <- as.character(cell_types[[name1]])
        cell_type2 <- as.character(cell_types[[name2]])

        #No need to calculate min distance when they're the same cell type
        if (name1 == name2) {
            local_result = data.frame(Target = name1, Nearest = name2, Mean = NA, Std.Dev = NA, Median = NA)
        } else {
            #vector to store all mins
            local_dist_min <- vector()
            #grab coordinates of all other cells in cell_type2 to calculate distances with
            all_celltype2_cord <- formatted_data[formatted_data$Phenotype == name2, c("Cell.X.Position", "Cell.Y.Position")]
            all_celltype1_cord <- formatted_data[formatted_data$Phenotype == name1, c("Cell.X.Position", "Cell.Y.Position")]

            #find all of closest points
            all_closest <- nn2(data = all_celltype2_cord, query = all_celltype1_cord, k = 1)

            local_dist_mins <- all_closest$nn.dists

            local_result <- data.frame(Target = name1, Nearest = name2, Mean = mean(local_dist_mins), Std.Dev = sd(local_dist_mins), Median = median(local_dist_mins))
        }
        result <- rbind(result, local_result)
    }

    # remove NAs e.g. for distance of cell against itself
    result <- result[complete.cases(result),]

    return(result)
}
