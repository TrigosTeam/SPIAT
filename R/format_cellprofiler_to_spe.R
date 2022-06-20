#' Format an image into a SpatialExperiment object
#'
#' @description Reads in spatial data in the form of cell coordinates, cell
#'   phenotypes (if available), and marker intensities and transforms to a
#'   SingleCellExperiment object. Users can input data in a general format
#'   (recommended), or data generated from inForm, HALO, Visum, CODEX,
#'   cellprofiler. The count assay stores the intensity level of every marker
#'   (rows) for every cell (columns). Cell phenotype, x and y coordinates and
#'   other cell properties are stored under colData. For the INFORM format, the
#'   cell properties needed to be included are Cell.Area, Nucleus.Area,
#'   Nucleus.Compactness, Nucleus.Axis.Ratio, and Cell.Axis.Ratio. For the HALO
#'   format, the cell properties to be included are Cell.Area, Nucleus.Area and
#'   Cytoplasm.Area. Note that if the data does not include these parameters, we
#'   recommend adding it to the output from inForm or HALO with NAs in columns,
#'   or alternatively using the ‘generic formatting’ option.
#'
#' @export
#' @param path String of the path location cellprofiler csv file.
#' @param markers String Vector containing the markers used for staining.
#' @param intensity_columns_interest String Vector with the names of the columns
#'   with the level of each marker. Column names must match the order of the
#'   'markers' parameter.
#' @return A SingleCellExperiment object is returned

format_cellprofiler_to_spe <- function(path = NULL, 
                                markers = NULL,  
                                intensity_columns_interest = NULL){
    
    Object.Id <- Cell.ID <- ObjectNumber <- NULL

    # Due to the implementation, rename DNA to DAPI before usage. 
    if (substr(path, nchar(path)-3+1, nchar(path)) != "csv"){
        stop("Please convert your file from a cpout to a csv format") }
    # read in the image
    image <- utils::read.csv(path)
    
    # CHECK 
    # if image contains all the columns specified and vectors of same length
    image_colnames <- colnames(image)
    if (!all(intensity_columns_interest %in% image_colnames)) {
        stop("One or more Intensity_columns_interest not found in image")}
        
    marker_count <- length(markers)
    intensity_col_count <- length(intensity_columns_interest)
    
    if (marker_count != intensity_col_count) {
        stop("The number of columns and markers do not match")
    }
        
    #extract intensities
    intensity_of_markers <- image[,c("ObjectNumber", intensity_columns_interest)]
    colnames(intensity_of_markers) <- c("ObjectNumber", markers)
    intensity_of_markers[intensity_of_markers == "#N/A"] <- NA
    intensity_of_markers <- remove_intensity_na(intensity_of_markers)
    intensity_of_markers <- apply(intensity_of_markers, 2, function(x){
        as.numeric(as.character(x))})
        
    # remove intensity NA rows from image
    image <- subset(image, ObjectNumber %in% intensity_of_markers[, "ObjectNumber"])
    intensity_of_markers <- intensity_of_markers[ , !(colnames(intensity_of_markers) == "ObjectNumber")]
    
    #grab relevant columns
    image <- image[,c("ObjectNumber", "Location_Center_X", "Location_Center_Y")]
    
    #rename Object.ID to Cell.ID
    colnames(image)[colnames(image) %in% c("ObjectNumber")] <- c("Cell.ID")

    #add "Cell_" in front of Cell.ID
    image$Cell.ID <- paste("Cell_", image$Cell.ID, sep="")
    
    colnames(image)[colnames(image) %in% c("Location_Center_X", "Location_Center_Y")] <- c("Cell.X.Position", "Cell.Y.Position")
    
    image$Phenotype <- NA
    image <- image[,c("Cell.ID", "Phenotype", "Cell.X.Position", "Cell.Y.Position")]
    
    #create the formatted_data with intensity levels
    formatted_data <- cbind(image, intensity_of_markers)
    
    #now create the SpE object...
    #grab the intensity level, markers and cell IDs
    assay_data <- formatted_data[,c(markers)]
    
    #transpose the matrix so every column is a cell and every row is a marker
    assay_data_matrix <- as.matrix(assay_data)
    colnames(assay_data_matrix) <- NULL
    rownames(assay_data_matrix) <- NULL
    assay_data_matrix_t <- t(assay_data_matrix)
    
    #Assign the phenotype, X and Y positions and cell property columns as the colData
    metadata_columns <- formatted_data[ ,c("Phenotype", "Cell.X.Position", 
                                           "Cell.Y.Position")]
    
    spe <- SpatialExperiment::SpatialExperiment(
        assay = assay_data_matrix_t,
        colData = metadata_columns,
        spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position"))
    
    rownames(spe) <- markers
    colnames(spe) <- formatted_data[,"Cell.ID"]

    return(spe)
}
