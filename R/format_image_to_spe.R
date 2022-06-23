#' Format an image into a SpatialExperiment object
#'
#' @description Reads in spatial data in the form of cell coordinates, cell
#'   phenotypes (if available), and marker intensities and transforms to a
#'   SpatialExperiment object. The assay stores the intensity level of every
#'   marker (rows) for every cell (columns). Cell phenotype is stored under
#'   `colData()`. Cell x and y coordinates are stored under `spatialCoords()`
#'   field. The function can read in data in general format (manually
#'   constructed input), or data from other platforms including inForm, HALO,
#'   CODEX and cellprofiler. Alternatively, users can use the specific function
#'   for each format.
#' @details Note for "cellprofiler" format, when specifying `markers`, please
#'   use "DAPI" to replace "DNA" due to implementation. The output data will
#'   include "DAPI" instead of "DNA".
#'
#'   The format of "Phenotype" column: For example, a cell positive for both
#'   "CD3" and "CD4" markers has the "CD3,CD4" **cell phenotype**. The phenotype
#'   has to be strictly formatted in such way where each positive marker has to
#'   be separated by a coma, with no space in between, and the order of the
#'   positive markers has to be the same as the order in the assay.
#'
#' @export
#' @param format String specifying the format of the data source. Default is
#'   "general" (RECOMMENDED), where the cell phenotypes, coordinates and marker
#'   intensities are imported manually by the user. Other formats include
#'   "inForm", "HALO", "cellprofiler" and "CODEX".
#' @param intensity_matrix (Optional) For "general" format. A matrix of marker
#'   intensities or gene expression where the column names are the Cell IDs, and
#'   the rownames the marker.
#' @param phenotypes (Optional) For "general" format. String Vector of cell
#'   phenotypes in the same order in which they appear in `intensity_matrix`. If
#'   no phenotypes available, then a vector of NAs can be used as input. Note
#'   that the combination of markers (e.g. CD3,CD4) needs to be used instead of
#'   the cell type name (e.g. helper T cells).
#' @param coord_x (Optional) For "general" format. Numeric Vector with the X
#'   coordinates of the cells. The cells must be in the same order as in the
#'   `intensity_matrix`.
#' @param coord_y (Optional) For "general" format. Numeric Vector with the Y
#'   coordinates of the cells. The cells must be in the same order as in the
#'   `intensity_matrix`.
#' @param path (Optional) For formats other than "general". String of the path
#'   location of the source file.
#' @param markers For formats other than "general". String Vector containing the
#'   markers used for staining. These must be in the same order as the marker
#'   columns in the input file, and must match the marker names used in the
#'   input file. One of the markers must be "DAPI".
#' @param locations (Optional) For "inForm" and "HALO". String Vector containing
#'   the locations of markers used for staining. Location can be either
#'   "Nucleus", "Cytoplasm" or "Membrane". This is used to select the Intensity
#'   column and can be used instead of `intensity_columns_interest`.
#' @param intensity_columns_interest (Optional) For "inForm" and "HALO", use if
#'   `locations` is not specified. For "cellprofiler", mandatory. Vector with
#'   the names of the columns with the level of each marker. Column names must
#'   match the order of the 'markers' parameter.
#' @param dye_columns_interest (Optional) For "HALO". Use if locations is not
#'   specified. Vector of names of the columns with the marker status (i.e.
#'   those indicating 1 or 0 for whether the cell is positive or negative for
#'   the marker). Column names must match the order of the 'markers' parameter.
#' @param path_to_codex_cell_phenotypes (Optional) For "CODEX".String of the
#'   path to the Cluster ID/Cell type file.
#' @seealso \code{\link{format_inform_to_spe}} \code{\link{format_halo_to_spe}}
#'   \code{\link{format_codex_to_spe}} \code{\link{format_cellprofiler_to_spe}}
#'
#' @return A SpatialExperiment object is returned
#' @examples
#' #Construct a marker intensity matrix (rows are markers, columns are cells)
#' intensity_matrix <- matrix(c(14.557, 0.169, 1.655, 0.054, 17.588, 0.229,
#' 1.188, 2.074, 21.262, 4.206,  5.924, 0.021), nrow = 4, ncol = 3)
#' # define marker names as rownames
#' rownames(intensity_matrix) <- c("DAPI", "CD3", "CD4", "AMACR")
#' # define cell IDs as colnames
#' colnames(intensity_matrix) <- c("Cell_1", "Cell_2", "Cell_3")
#' # Construct a dummy metadata (phenotypes, x/y coordinates)
#' # the order of the elements in these vectors correspond to the cell order
#' # in `intensity matrix`
#' phenotypes <- c("OTHER",  "AMACR", "CD3,CD4")
#' coord_x <- c(82, 171, 184)
#' coord_y <- c(30, 22, 38)
#'
#' formatted_image <- format_image_to_spe(intensity_matrix=intensity_matrix,
#' phenotypes = phenotypes, coord_x = coord_x,coord_y = coord_y)

format_image_to_spe <- function(format = "general", intensity_matrix = NULL, 
                                phenotypes = NULL, coord_x = NULL, 
                                coord_y = NULL, path = NULL, markers = NULL, 
                                locations = NULL,
                                intensity_columns_interest = NULL, 
                                dye_columns_interest = NULL,
                                path_to_codex_cell_phenotypes = NULL){
    if (format == "general"){
        intensity_matrix <- remove_intensity_na(intensity_matrix)
        
        metadata_columns <- data.frame(Phenotype = phenotypes,
                                       Cell.X.Position = coord_x,
                                       Cell.Y.Position = coord_y)
        spe <- SpatialExperiment::SpatialExperiment(
            assay = intensity_matrix,
            colData = metadata_columns,
            spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position"))
    } else if (format == "inForm") {
       spe <- format_inform_to_spe(path = path, markers = markers,
                                   locations = locations, 
                                   intensity_columns_interest = 
                                       intensity_columns_interest)
    }else if (format == "HALO"){
       spe <- format_halo_to_spe(path = path, markers = markers,
                                 locations = locations, 
                                 dye_columns_interest = dye_columns_interest,
                                 intensity_columns_interest = 
                                     intensity_columns_interest)
    } else if(format == "CODEX"){
       spe <- format_codex_to_spe(path = path, markers = markers, 
                 path_to_codex_cell_phenotypes = path_to_codex_cell_phenotypes)
    }else if(format == "cellprofiler") {
       spe <- format_cellprofiler_to_spe(path = path, markers = markers, 
                 intensity_columns_interest = intensity_columns_interest)
    } else { methods::show("Please entre a valid format!" )}
    return(spe)
}
