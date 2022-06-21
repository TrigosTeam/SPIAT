#' Format a CODEX image into a SpatialExperiment object
#'
#' @description Reads in spatial data in the form of cell coordinates, cell
#'   phenotypes (if available), and marker intensities and transforms to a
#'   `SpatialExperiment` object. The assay stores the intensity level of every marker
#'   (rows) for every cell (columns). Cell phenotype is stored under colData. 
#'   Cell x and y coordinates are stored under `spatialCoords()` field. 
#' @export
#' @param path String of the path location of CODEX csv file. 
#' @param markers String Vector containing the markers used for staining.
#' @param path_to_codex_cell_phenotypes String of the path to
#'   the Cluster ID/Cell type file.
#' @return A SpatialExperiment object is returned
#' @examples 
#' path <- system.file("extdata", "tiny_codex.csv.gz", package = "SPIAT")
#' path_to_codex_cell_phenotypes <- system.file("extdata", 
#' "tiny_codex_phenotypes.txt.gz", package = "SPIAT")
#' markers <- c("CD45", "Ly6C", "CD27", "CD5", "CD79b")
#' formatted_codex <- format_codex_to_spe(path = path, markers = markers,
#' path_to_codex_cell_phenotypes = path_to_codex_cell_phenotypes)

format_codex_to_spe <- function(path = NULL, markers,
                                path_to_codex_cell_phenotypes = NULL){
    phenotypes <- utils::read.delim(path_to_codex_cell_phenotypes,header=FALSE)
    colnames(phenotypes) <- c("Cluster_ID", "Cell_type")
    
    data <- utils::read.delim(path, sep=",")
    data$Cell_type <- phenotypes$Cell_type[match(data$Imaging.phenotype.cluster.ID,
                                                 phenotypes$Cluster_ID)]
    data$Imaging.phenotype.cluster.ID <- NULL
    data$niche.cluster.ID <- NULL  
    data$sample_Xtile_Ytile <- NULL
    data$Z.Z <- NULL
    
    coordinates <- data[,c("X.X", "Y.Y")]
    data$X.X <- NULL
    data$Y.Y <- NULL
    phenotype <- data$Cell_type
    data$Cell_type <- NULL
    rownames(data) <- paste("Cell", rownames(data), sep="_")
    data <- t(data)

    metadata_columns <- data.frame(Phenotype = phenotype,
                                   Cell.X.Position = coordinates$X.X,
                                   Cell.Y.Position = coordinates$Y.Y)                     
    
    spe <- SpatialExperiment::SpatialExperiment(
        assay = data,
        colData = metadata_columns,
        spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position"))
    
    rownames(spe) <- markers
    
    return(spe)
}
