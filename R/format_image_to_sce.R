#' format_image_to_sce
#'
#' @description Formats an INFORM or HALO image into a singlecellexperiment class
#' where the count assay stores the expression level of every marker (rows) for
#' every cell (columns), and cell phenotype, x and y coordinates are stored
#' under colData
#'
#' @export
#' @param format String defining the software use for cell segmentation. Options: "INFORM" or "HALO"
#' @param image String of the path location of either HALO csv file or INFORM textfile
#' @param markers Vector containing the markers used for staining
#' @param dye_columns_interest (Only for HALO formats) Vector of names of the columns with the marker status 
#' (i.e. those indicating 1 or 0 for whether the cell is positive or negative for the marker). 
#' Column names must match the order of the 'markers' parameter.
#' @param intensity_columns_interest Vector with the names of the columns with the level of each marker.
#' Column names must match the order of the 'markers' parameter
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom utils read.csv read.delim

format_image_to_sce <- function(format = "INFORM", image, markers, dye_columns_interest = NULL, intensity_columns_interest) {

    #replace the spaces and non-alphanumeric characters as a '.' for column selection
    intensity_columns_interest <- gsub("[^[:alnum:]]", ".", intensity_columns_interest)
    dye_columns_interest <- gsub("[^[:alnum:]]", ".", dye_columns_interest)
    
    #process the data based on data format
    if (format == "HALO"){
        #following is from format_HALO_new

        #read in the image
        image <- read.csv(image)
        
        #CHECK - if image contains all the columns specified and vectors of same length
        image_colnames <- colnames(image)
        if (!all(intensity_columns_interest %in% image_colnames)) {
            stop("One or more Intensity_columns_interest not found in image")
        }
        if (!all(dye_columns_interest %in% image_colnames)) {
            stop("One or more dye_columns_interest not found in image")
        }
        marker_count <- length(markers)
        intensity_col_count <- length(intensity_columns_interest)
        dye_col_count <- length(dye_columns_interest)
        if (marker_count != intensity_col_count || marker_count != dye_col_count || intensity_col_count != dye_col_count) {
            stop("The number of dyes, columns and markers do not match")
        }

        #First remove all non-DAPI cells
        idx <- which(markers == "DAPI")
        
        #CHECK - if DAPI is in the dyes
        if (length(idx) == 0) {
            stop("Please include DAPI in the markers")
        }
        
        DAPI_col_name <- dye_columns_interest[idx]
        DAPI_non_zero_rows <- which(image[,DAPI_col_name] != 0)
        image <- image[DAPI_non_zero_rows,]
        
        #extract intensities
        intensity_of_markers <- image[,intensity_columns_interest]
        colnames(intensity_of_markers) <- markers
        intensity_of_markers[intensity_of_markers == "#N/A"] <- NA
        intensity_of_markers <- apply(intensity_of_markers, 2, function(x){
            as.numeric(as.character(x))
        })

        #get the expression status columns
        expression_status_cols <- image[,dye_columns_interest]
        colnames(expression_status_cols) <- markers

        #grab relavant columns
        image <- image[,c("Object.Id", "XMin", "XMax", "YMin", "YMax")]

        #rename Object.ID to Cell.ID
        colnames(image)[1] <- "Cell.ID"

        #add "Cell_" in front of Cell.ID
        image$Cell.ID <- paste("Cell_", image$Cell.ID, sep="")

        #add averaged X and Y position
        image$Cell.X.Position <- (image$XMin + image$XMax)/2
        image$Cell.Y.Position <- (image$YMin + image$YMax)/2

        #start reading in the Phenotypes of every cell
        expression_status_cols$Phenotype <- ""
        for (marker in markers) {
            if (marker == "DAPI") {
                phenotype <- "OTHER,"
            } else {
                phenotype <- paste(marker, ",", sep = "")
            }

            #get the row idx of the cells that express the specific marker, and paste the phenotype
            rows_true_exp <- which(expression_status_cols[,marker] != 0)
            if (length(rows_true_exp) != 0) {
                expression_status_cols[rows_true_exp,]$Phenotype <- paste(expression_status_cols[rows_true_exp,]$Phenotype, phenotype, sep="")
            }
        }

        #now clean the phenotype column
        if (nrow(expression_status_cols[expression_status_cols$Phenotype == "OTHER,", ]) != 0) {
            expression_status_cols[expression_status_cols$Phenotype == "OTHER,", ]$Phenotype <- "OTHER"
        }
        expression_status_cols$Phenotype <- gsub("OTHER,", "", expression_status_cols$Phenotype)
        expression_status_cols$Phenotype <- gsub(",OTHER", "", expression_status_cols$Phenotype)
        expression_status_cols$Phenotype <- gsub(",$", "", expression_status_cols$Phenotype)

        #grab the phenotype column and cbind to image
        phenotype_column <- data.frame(expression_status_cols$Phenotype)
        colnames(phenotype_column) <- "Phenotype"
        image <- cbind(image, phenotype_column)

        image$Phenotype <- as.character(image$Phenotype)

        image <- image[,c("Cell.ID", "Phenotype", "Cell.X.Position", "Cell.Y.Position")]

    } else if (format == "INFORM"){

        ###following codes are from format_INFORM

        #open the image file
        image <- read.delim(image)
        #remove all rows with empty phenotype/no markers
        image <- image[image$Phenotype != "",]
        image <- image[!is.na(image$Phenotype), ]
        
        #CHECK - if image contains all the columns specified and vectors of same length
        image_colnames <- colnames(image)
        if (!all(intensity_columns_interest %in% image_colnames)) {
            stop("One or more Intensity_columns_interest not found in image")
        }
        marker_count <- length(markers)
        intensity_col_count <- length(intensity_columns_interest)
        if (marker_count != intensity_col_count) {
            stop("The number of dyes and columns does not match")
        }

        ###added: extract intensities
        intensity_of_markers <- image[,intensity_columns_interest]
        colnames(intensity_of_markers) <- markers
        intensity_of_markers[intensity_of_markers == "#N/A"] <- NA
        intensity_of_markers <- apply(intensity_of_markers, 2, function(x){
            as.numeric(as.character(x))
        })

        #extract the columns of interest and discard the rest
        image <- image[ ,c("Cell.ID", "Phenotype", "Cell.X.Position",
                           "Cell.Y.Position")]

        #add 'cell_' to the start of the objectId
        image$Cell.ID <- paste0("Cell_", image$Cell.ID)

        #standardize PDL-1 into PDL1
        image$Phenotype <- gsub("PDL-1", "PDL1", image$Phenotype, fixed=TRUE)

        #reformat the phenotype into "marker1, marker2..."
        for (marker in markers) {
            marker_positive <- paste0(marker, "+", sep="")
            marker_negative <- paste0(marker, "-", sep="")

            marker_replacement <- paste0(marker,",", sep="")

            image$Phenotype <- gsub(marker_positive, marker_replacement, image$Phenotype, fixed=TRUE)
            image$Phenotype <- gsub(marker_negative, "", image$Phenotype, fixed=TRUE)

        }

        #remove the comma at the end
        image$Phenotype <- gsub(",$","", image$Phenotype)

    } else{
        stop("Please enter a valid format: INFORM/HALO")
    }

    #create the formatted_data with intensity levels
    formatted_data <- cbind(image, intensity_of_markers)

    #now create the SCE object...
    #grab the expression level, markers and cell IDs
    assay_data <- formatted_data[,markers]
    assay_rownames <- markers
    assay_colnames <- formatted_data[,"Cell.ID"]

    #transpose the matrix so every column is a cell and every row is a marker
    assay_data_matrix <- as.matrix(assay_data)
    colnames(assay_data_matrix) <- NULL
    rownames(assay_data_matrix) <- NULL
    assay_data_matrix_t <- t(assay_data_matrix)

    sce <- SingleCellExperiment(assays = list(counts = assay_data_matrix_t))

    rownames(sce) <- assay_rownames
    colnames(sce) <- assay_colnames

    #Assign the phenotype, X and Y positions as the colData
    coldata_phenotype <- formatted_data[,"Phenotype"]
    coldata_Xpos <- formatted_data[,"Cell.X.Position"]
    coldata_Ypos <- formatted_data[,"Cell.Y.Position"]
    SummarizedExperiment::colData(sce)$Phenotype <- coldata_phenotype
    SummarizedExperiment::colData(sce)$Cell.X.Position <- coldata_Xpos
    SummarizedExperiment::colData(sce)$Cell.Y.Position <- coldata_Ypos

    return(sce)
}


