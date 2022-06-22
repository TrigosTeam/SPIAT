#' Format an inForm image into a SpatialExperiment object
#'
#' @description Reads in inForm data in the form of cell coordinates, cell
#'   phenotypes (if available), and marker intensities and transforms to a
#'   SpatialExperiment object.  The assay stores the intensity level of
#'   every marker (rows) for every cell (columns). Cell phenotype, x and y
#'   coordinates and other cell properties are stored under colData. The cell
#'   properties to include are Cell.Area, Nucleus.Area, Nucleus.Compactness,
#'   Nucleus.Axis.Ratio, and Cell.Axis.Ratio.  Note that if the data does not
#'   include these parameters, we recommend adding it to the output from inForm
#'   with NAs in columns.
#'
#' @export
#' @param path String of the path location of inForm text file.
#' @param markers String Vector containing the markers used for staining.
#' @param locations (Optional) String Vector containing the locations of markers
#'   used for staining. Location can be either "Nucleus", "Cytoplasm" or
#'   "Membrane". This is used to select the Intensity column and can be used
#'   instead of `intensity_columns_interest`.
#' @param intensity_columns_interest (Optional) Use if `locations` is not
#'   specified. Vector with the names of the columns with the level of each
#'   marker. Column names must match the order of the 'markers' parameter.
#' @return A SpatialExperiment object is returned
#' @examples
#' raw_inform_data<-system.file("extdata","tiny_inform.txt.gz",package="SPIAT")
#' markers <- c("DAPI", "CD3", "PD-L1", "CD4", "CD8", "AMACR")
#' locations <- c("Nucleus", "Cytoplasm", "Membrane", "Cytoplasm", "Cytoplasm", 
#' "Cytoplasm")
#' formatted_inForm <- format_inform_to_spe(path=raw_inform_data, 
#' markers=markers, locations=locations)

format_inform_to_spe <- function(path, markers, locations=NULL,  
                                intensity_columns_interest=NULL){
    
    Object.Id <- Cell.ID <- ObjectNumber <- NULL

    #read in the image file
    image <- vroom::vroom(path)
    #if Phenotype column exists,
    #remove all rows with empty phenotype/no markers
    if (is.null(image$Phenotype)){
        image <- image[image$Phenotype != "",]
        image <- image[!is.na(image$Phenotype), ]
    }
    # rename "Cell ID" to "Cell.ID"
    colnames(image)[which(names(image) == "Cell ID")] <- "Cell.ID"
    
    if (!is.null(locations)) {
        # add the location of interest to each marker
        names_to_match <- paste(locations, markers, sep=" ")
        
        # get all the mean intensity column names in the file
        intensity_col_all <- 
            colnames(image)[
                grepl("Mean \\(Normalized Counts, Total Weighting\\)",
                                  colnames(image))]
        
        # get the intensity column name for each marker
        intensity_columns_interest <- character(length(markers))
        i <- 1
        for (name in names_to_match) {
            intensity_columns_interest[i] <- 
                intensity_col_all[grepl(name, intensity_col_all)]
            i <- i + 1
        }
    } else {
        #CHECK
        image_colnames <- colnames(image)
        if (!all(intensity_columns_interest %in% image_colnames)) {
            stop("One or more Intensity_columns_interest not found in image")
        }
        marker_count <- length(markers)
        intensity_col_count <- length(intensity_columns_interest)
        if (marker_count != intensity_col_count) {
            stop("The number of dyes and columns does not match")
        }
    }
        
    ###added: extract intensities
    intensity_of_markers <- image[,c("Cell.ID", intensity_columns_interest)]
    colnames(intensity_of_markers) <- c("Cell.ID", markers)
    intensity_of_markers[intensity_of_markers == "#N/A"] <- NA
    intensity_of_markers <- remove_intensity_na(intensity_of_markers)
    intensity_of_markers <- apply(intensity_of_markers, 2, function(x){
        as.numeric(as.character(x))
    })
    
    # remove intensity NA rows from image
    image <- subset(image, `Cell.ID` %in% intensity_of_markers[, "Cell.ID"])
    intensity_of_markers <- 
        intensity_of_markers[ , !(colnames(intensity_of_markers) == "Cell.ID")]
    
    #extract the columns of interest and discard the rest
    colnames(image) <- make.names(colnames(image))
    
    # shorten column names of some properties 
    names(image)[names(image)=="Entire.Cell.Area..pixels."] <- "Cell.Area"
    names(image)[names(image)=="Nucleus.Area..pixels."] <- "Nucleus.Area"
    names(image)[names(image)=="Entire.Cell.Axis.Ratio"] <- "Cell.Axis.Ratio"
    
    cell_properties_cols <- c("Cell.Area","Nucleus.Area",
                              "Nucleus.Compactness","Nucleus.Axis.Ratio",
                              "Cell.Axis.Ratio")
    image <- image[ ,c("Cell.ID", "Phenotype", "Cell.X.Position",
                       "Cell.Y.Position", cell_properties_cols)]
    
    #add 'cell_' to the start of the objectId
    image$Cell.ID <- paste0("Cell_", image$Cell.ID)
        
    #reformat the phenotype into "marker1, marker2..."
    for (marker in markers) {
        marker_positive <- paste0(marker, "+", sep="")
        marker_negative <- paste0(marker, "-", sep="")
        
        marker_replacement <- paste0(marker,",", sep="")
        
        image$Phenotype <- gsub(marker_positive, marker_replacement, 
                                image$Phenotype, fixed=TRUE)
        image$Phenotype <- gsub(marker_negative, "", image$Phenotype, 
                                fixed=TRUE)
    }
    
    #remove the comma at the end
    image$Phenotype <- gsub(",$","", image$Phenotype)
    
    #create the formatted_data with intensity levels
    formatted_data <- cbind(image, intensity_of_markers)
    
    #now create the spe object...
    #grab the intensity level, markers and cell IDs
    assay_data <- formatted_data[,c(markers)]
    
    #transpose the matrix so every column is a cell and every row is a marker
    assay_data_matrix <- as.matrix(assay_data)
    colnames(assay_data_matrix) <- NULL
    rownames(assay_data_matrix) <- NULL
    assay_data_matrix_t <- t(assay_data_matrix)
    
    #Assign the phenotype, X and Y positions and cell property columns
    metadata_columns <- formatted_data[,c("Phenotype", "Cell.X.Position", 
                                          "Cell.Y.Position", 
                                          cell_properties_cols)]
    
    spe <- SpatialExperiment::SpatialExperiment(
        assay = assay_data_matrix_t,
        colData = metadata_columns,
        spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position"))
        
    rownames(spe) <- markers
    colnames(spe) <- formatted_data[,"Cell.ID"]
    
    return(spe)
}
