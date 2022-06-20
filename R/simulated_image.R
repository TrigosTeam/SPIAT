#' SPE object of a formatted image (simulated by `spaSim` package)
#'
#' A dataset that contains a formatted spe object with cell ids and phenotypes
#' in `colData()` and marker intensities in `assays()`.
#'
#' @format An SpatialExperiment object. Assay contains 5 rows (markers) and 4951
#'   columns (cells); colData contains 4951 rows (cells) and 3 columns.
#' @seealso \code{\link{defined_image}} \code{\link{image_no_markers}}
"simulated_image"
