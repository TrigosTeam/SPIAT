#' Format a HALO image into a SpatialExperiment object
#'
#' @description Reads in HALO data in the form of cell coordinates, cell
#'   phenotypes (if available), and marker intensities and transforms to a
#'   `SpatialExperiment` object. The assay stores the intensity level of every
#'   marker (rows) for every cell (columns). Cell x and y coordinates are stored
#'   under `spatialCoords()`. Cell phenotype and other cell properties are
#'   stored under colData. The cell properties to be included are Cell.Area,
#'   Nucleus.Area and Cytoplasm.Area. Note that if the data does not include
#'   these parameters, we recommend adding it to the output from HALO with NAs
#'   in columns.
#'
#' @export
#' @param path String of the path location of HALO text file.
#' @param markers String Vector containing the markers used for staining.
#' @param locations (Optional) Vector containing the locations of markers used for
#'   staining. Location can be either "Nucleus", "Cytoplasm" or "Membrane". This
#'   is used to select the Intensity column and can be used instead of
#'   `intensity_columns_interest`.
#' @param dye_columns_interest (Optional) Use if locations is not
#'   specified. Vector of names of the columns with the marker status (i.e.
#'   those indicating 1 or 0 for whether the cell is positive or negative for
#'   the marker). Column names must match the order of the 'markers' parameter.
#' @param intensity_columns_interest (Optional) Use if locations is not
#'   specified. Vector with the names of the columns with the level of each
#'   marker. Column names must match the order of the 'markers' parameter.
#' @return A SpatialExperiment object is returned
#' @examples
#' raw_halo_data <- system.file("extdata", "tiny_halo.csv.gz", package="SPIAT")
#' markers <- c("DAPI", "CD3", "PDL-1", "CD4", "CD8", "AMACR")
#' intensity_columns_interest <- c("Dye 1 Nucleus Intensity",
#' "Dye 2 Cytoplasm Intensity","Dye 3 Membrane Intensity",
#' "Dye 4 Cytoplasm Intensity", "Dye 5 Cytoplasm Intensity",
#' "Dye 6 Cytoplasm Intensity")
#' dye_columns_interest <-c("Dye 1 Positive Nucleus","Dye 2 Positive Cytoplasm",
#' "Dye 3 Positive Membrane", "Dye 4 Positive Cytoplasm",
#' "Dye 5 Positive Cytoplasm", "Dye 6 Positive Cytoplasm")
#' formatted_HALO <- format_halo_to_spe(path=raw_halo_data,markers=markers, 
#' intensity_columns_interest=intensity_columns_interest,
#' dye_columns_interest=dye_columns_interest)

format_halo_to_spe <- function(path = NULL, 
                                markers = NULL, 
                                locations = NULL, 
                                dye_columns_interest = NULL, 
                                intensity_columns_interest = NULL){
    
    Object.Id <- Cell.ID <- ObjectNumber <- NULL
    #read in the image using vroom for super fast import
    image <- vroom::vroom(path)
    
    #remove spaces from column names
    colnames(image) <- make.names(colnames(image))
    
    # if locations is specified, use location plus marker name to get the intensity and dye columns
    if (!is.null(locations)) {
        intensity_columns_interest <- character(length(markers))
        dye_columns_interest <- character(length(markers))
        i <- 1
            
        for (loc in locations) {
            if (loc == "Nucleus") {
                intensity_columns_interest[i] <- 
                    paste0("Dye.", i, ".Nucleus.Intensity")
                dye_columns_interest[i] <- 
                    paste0("Dye.", i, ".Positive.Nucleus")
            } else if (loc == "Cytoplasm") {
                intensity_columns_interest[i] <- 
                    paste0("Dye.", i, ".Cytoplasm.Intensity")
                dye_columns_interest[i] <- 
                    paste0("Dye.", i, ".Positive.Cytoplasm")
            } else if (loc == "Membrane") {
                intensity_columns_interest[i] <- 
                    paste0("Dye.", i, ".Membrane.Intensity")
                dye_columns_interest[i] <- 
                    paste0("Dye.", i, ".Positive.Membrane")
            } else {
                stop('Location incorrectly specified. Must be either "Nucleus", 
                     "Cytoplasm" or "Membrane"')
            }
            i <- i + 1}
            
    # if locations not used, get the intensity and dye columns from their specified names
    } else {
        #replace the spaces and non-alphanumeric characters as a '.' for column selection
        intensity_columns_interest <- gsub("[^[:alnum:]]", ".", 
                                           intensity_columns_interest)
        dye_columns_interest <- gsub("[^[:alnum:]]", ".", dye_columns_interest)
        
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
        if (marker_count != intensity_col_count || marker_count != 
            dye_col_count || intensity_col_count != dye_col_count) {
            stop("The number of dyes, columns and markers do not match")
        }
    }
        
    #First remove all non-DAPI cells
    idx <- which(markers == "DAPI")
    #CHECK - if DAPI is in the dyes
    if (length(idx) == 0) {
        stop("Please include DAPI in the markers")}
    
    DAPI_col_name <- dye_columns_interest[idx]
    DAPI_non_zero_rows <- which(image[,DAPI_col_name] != 0)
    image <- image[DAPI_non_zero_rows,]
    
    #extract intensities
    intensity_of_markers <- image[, c("Object.Id", intensity_columns_interest)]
    colnames(intensity_of_markers) <- c("Object.Id", markers)
    intensity_of_markers[intensity_of_markers == "#N/A"] <- NA
    intensity_of_markers <- remove_intensity_na(intensity_of_markers)
    intensity_of_markers <- apply(intensity_of_markers, 2, function(x){
        as.numeric(as.character(x))})
        
    # remove intensity NA rows from image
    image <- subset(image, Object.Id %in% intensity_of_markers[, "Object.Id"])
    intensity_of_markers <- 
        intensity_of_markers[ , !(colnames(intensity_of_markers) == "Object.Id")]
    
    #get the intensity status columns
    intensity_status_cols <- image[,dye_columns_interest]
    colnames(intensity_status_cols) <- markers
    
    #grab relevant columns
    cell_properties_cols <- 
        colnames(image)[grep("Nucleus.Area|Cytoplasm.Area|Cell.Area", 
                             colnames(image))]
    image <- image[,c("Object.Id", "XMin", "XMax", "YMin", "YMax", 
                      cell_properties_cols)]
    
    #rename Object.ID to Cell.ID
    colnames(image)[1] <- "Cell.ID"
    
    #add "Cell_" in front of Cell.ID
    image$Cell.ID <- paste("Cell_", image$Cell.ID, sep="")
    
    #add averaged X and Y position
    image$Cell.X.Position <- (image$XMin + image$XMax)/2
    image$Cell.Y.Position <- (image$YMin + image$YMax)/2
    
    #start reading in the Phenotypes of every cell
    intensity_status_cols$Phenotype <- ""
    for (marker in markers) {
        if (marker == "DAPI") {
            phenotype <- "OTHER,"} else {
            phenotype <- paste(marker, ",", sep = "")}
        
        #get the row idx of the cells that express the specific marker, and paste the phenotype
        rows_true_exp <- which(intensity_status_cols[,marker] != 0)
        if (length(rows_true_exp) != 0) {
            intensity_status_cols[rows_true_exp,]$Phenotype <- 
                paste(intensity_status_cols[rows_true_exp,]$Phenotype, 
                      phenotype, sep="")
        }}
        
    #now clean the phenotype column
    if (nrow(intensity_status_cols[intensity_status_cols$Phenotype == "OTHER,", ]) != 0) {
        intensity_status_cols[intensity_status_cols$Phenotype == "OTHER,", ]$Phenotype <- "OTHER"
    }
    intensity_status_cols$Phenotype <- gsub("OTHER,", "", 
                                            intensity_status_cols$Phenotype)
    intensity_status_cols$Phenotype <- gsub(",OTHER", "", 
                                            intensity_status_cols$Phenotype)
    intensity_status_cols$Phenotype <- gsub(",$", "", 
                                            intensity_status_cols$Phenotype)
    
    #grab the phenotype column and cbind to image
    phenotype_column <- data.frame(intensity_status_cols$Phenotype)
    colnames(phenotype_column) <- "Phenotype"
    
    image <- cbind(image, phenotype_column)
    image$Phenotype <- as.character(image$Phenotype)
    image <- image[,c("Cell.ID", "Phenotype", "Cell.X.Position", 
                      "Cell.Y.Position", cell_properties_cols)]
    
    #create the formatted_data with intensity levels
    formatted_data <- cbind(image, intensity_of_markers)
    
    #now create the SPE object...
    #grab the intensity level, markers and cell IDs
    assay_data <- formatted_data[,c(markers)]
    
    #transpose the matrix so every column is a cell and every row is a marker
    assay_data_matrix <- as.matrix(assay_data)
    colnames(assay_data_matrix) <- NULL
    rownames(assay_data_matrix) <- NULL
    assay_data_matrix_t <- t(assay_data_matrix)
    
    #Assign the phenotype, X and Y positions and cell property columns as the colData
    metadata_columns <- formatted_data[,c("Phenotype", "Cell.X.Position", 
                                    "Cell.Y.Position", cell_properties_cols)]
    spe <- SpatialExperiment::SpatialExperiment(
        assay = assay_data_matrix_t,
        colData = metadata_columns,
        spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position"))
    
    rownames(spe) <- markers
    colnames(spe) <- formatted_data[,"Cell.ID"]
    
    return(spe)
}
