#' format_image_to_sce
#'
#' @description Formats an INFORM or HALO image into a singlecellexperiment class
#' where the count assay stores the intensity level of every marker (rows) for
#' every cell (columns). Cell phenotype, x and y coordinates and other properties are stored under colData. If InForm format, the cell properties 
#' are Cell.Area, Nucleus.Area, Nucleus.Compactness, Nucleus.Axis.Ratio, and Cell.Axis.Ratio. If HALO format, the cell properties are Cell.Area, 
#' Nucleus.Area, Cytoplasm.Area and Membrane.Perimeter.
#' 
#' @export
#' @param format String defining the software use for cell segmentation. Options: "INFORM" or "HALO"
#' @param image String of the path location of either HALO csv file or INFORM textfile
#' @param markers Vector containing the markers used for staining. Must match the order of the 'markers' parameter
#' @param locations Vector containing the locations of markers used for staining. Location can be either "Nucleus", "Cytoplasm" or "Membrane". 
#' This is used to select the Intensity column and can be used instead of intensity_columns_interest. 
#' @param dye_columns_interest (Only for HALO formats) Use if locations is not specified. Vector of names of the columns with the marker status 
#' (i.e. those indicating 1 or 0 for whether the cell is positive or negative for the marker). 
#' Column names must match the order of the 'markers' parameter.
#' @param intensity_columns_interest Use if locations is not specified. Vector with the names of the columns with the level of each marker.
#' Column names must match the order of the 'markers' parameter
#' @importFrom SingleCellExperiment SingleCellExperiment	
#' @importFrom SummarizedExperiment colData
#' @importFrom utils read.csv read.delim

format_image_to_sce <- function(format = "INFORM", image, markers, locations = NULL, dye_columns_interest = NULL, intensity_columns_interest = NULL) {
  
  #process the data based on data format
  if (format == "HALO"){
    #following is from format_HALO_new
    
    #read in the image
    image <- read.csv(image)
    
    # if locations is specified, use location plus marker name to get the intensity and dye columns
    if (!is.null(locations)) {
      
      intensity_columns_interest <- character(length(markers))
      dye_columns_interest <- character(length(markers))
      i <- 1
      
      for (loc in locations) {
        if (loc == "Nucleus") {
          intensity_columns_interest[i] <- paste0("Dye.", i, ".Nucleus.Intensity")
          dye_columns_interest[i] <- paste0("Dye.", i, ".Positive.Nucleus")
        } else if (loc == "Cytoplasm") {
          intensity_columns_interest[i] <- paste0("Dye.", i, ".Cytoplasm.Intensity")
          dye_columns_interest[i] <- paste0("Dye.", i, ".Positive.Cytoplasm")
        } else if (loc == "Membrane") {
          intensity_columns_interest[i] <- paste0("Dye.", i, ".Membrane.Intensity")
          dye_columns_interest[i] <- paste0("Dye.", i, ".Positive.Membrane")
        } else {
          stop('Location incorrectly specified. Must be either "Nucleus", "Cytoplasm" or "Membrane"')
        }
        
        i <- i + 1
      }
      
      # if locations not used, get the intensity and dye columns from their specified names
    } else {
      
      #replace the spaces and non-alphanumeric characters as a '.' for column selection
      intensity_columns_interest <- gsub("[^[:alnum:]]", ".", intensity_columns_interest)
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
    
    #get the intensity status columns
    intensity_status_cols <- image[,dye_columns_interest]
    colnames(intensity_status_cols) <- markers
    
    #grab relevant columns
    cell_properties_cols <- c("Nucleus.Area", "Cytoplasm.Area", "Membrane.Perimeter", "Cell.Area")
    image <- image[,c("Object.Id", "XMin", "XMax", "YMin", "YMax", cell_properties_cols)]
    
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
        phenotype <- "OTHER,"
      } else {
        phenotype <- paste(marker, ",", sep = "")
      }
      
      #get the row idx of the cells that express the specific marker, and paste the phenotype
      rows_true_exp <- which(intensity_status_cols[,marker] != 0)
      if (length(rows_true_exp) != 0) {
        intensity_status_cols[rows_true_exp,]$Phenotype <- paste(intensity_status_cols[rows_true_exp,]$Phenotype, phenotype, sep="")
      }
    }
    
    #now clean the phenotype column
    if (nrow(intensity_status_cols[intensity_status_cols$Phenotype == "OTHER,", ]) != 0) {
      intensity_status_cols[intensity_status_cols$Phenotype == "OTHER,", ]$Phenotype <- "OTHER"
    }
    intensity_status_cols$Phenotype <- gsub("OTHER,", "", intensity_status_cols$Phenotype)
    intensity_status_cols$Phenotype <- gsub(",OTHER", "", intensity_status_cols$Phenotype)
    intensity_status_cols$Phenotype <- gsub(",$", "", intensity_status_cols$Phenotype)
    
    #grab the phenotype column and cbind to image
    phenotype_column <- data.frame(intensity_status_cols$Phenotype)
    colnames(phenotype_column) <- "Phenotype"
    image <- cbind(image, phenotype_column)
    
    image$Phenotype <- as.character(image$Phenotype)
    
    image <- image[,c("Cell.ID", "Phenotype", "Cell.X.Position", "Cell.Y.Position", cell_properties_cols)]
    
  } else if (format == "INFORM"){
    
    ###following codes are from format_INFORM
    
    #open the image file, keep column names as is for matching to markers
    image <- read.delim(image, check.names=FALSE)
    #remove all rows with empty phenotype/no markers
    image <- image[image$Phenotype != "",]
    image <- image[!is.na(image$Phenotype), ]
    
    if (!is.null(locations)) {
    
      # add the location of interest to each marker
      names_to_match <- paste(locations, markers, sep=" ")
      
      # get all the mean intensity column names in the file
      intensity_col_all <- colnames(image)[grepl("Mean \\(Normalized Counts, Total Weighting\\)", colnames(image))]
      
      # get the intensity column name for each marker
      intensity_columns_interest <- character(length(markers))
      i <- 1
      
      for (name in names_to_match) {
        intensity_columns_interest[i] <- intensity_col_all[grepl(name, intensity_col_all)]
        i <- i + 1
      }
    
      } else {
    
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
    }
    
    ###added: extract intensities
    intensity_of_markers <- image[,intensity_columns_interest]
    colnames(intensity_of_markers) <- markers
    intensity_of_markers[intensity_of_markers == "#N/A"] <- NA
    intensity_of_markers <- apply(intensity_of_markers, 2, function(x){
      as.numeric(as.character(x))
    })
    
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
    
    #standardize PDL-1 into PDL1
    #image$Phenotype <- gsub("PDL-1", "PDL1", image$Phenotype, fixed=TRUE)
    
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
  #grab the intensity level, markers and cell IDs
  assay_data <- formatted_data[,c(markers)]
  assay_rownames <- c(markers)
  assay_colnames <- formatted_data[,"Cell.ID"]
  
  #transpose the matrix so every column is a cell and every row is a marker
  assay_data_matrix <- as.matrix(assay_data)
  colnames(assay_data_matrix) <- NULL
  rownames(assay_data_matrix) <- NULL
  assay_data_matrix_t <- t(assay_data_matrix)
  
  sce <- SingleCellExperiment(assays = list(counts = assay_data_matrix_t))
  
  rownames(sce) <- assay_rownames
  colnames(sce) <- assay_colnames
  
  #Assign the phenotype, X and Y positions and cell property columns as the colData
  metadata_columns <- formatted_data[ ,c("Phenotype", "Cell.X.Position", "Cell.Y.Position", cell_properties_cols)]
  SummarizedExperiment::colData(sce) <- cbind(SummarizedExperiment::colData(sce), metadata_columns)
  
  return(sce)
}
