#' SCE object of a formatted image (simulated by `spaSim` package)
#'
#' A dataset that contains a formatted sce object with cell ids, cell positions,
#' and phenotypes in metadata and marker intensities in assays.
#'
#' @format An sce object. Assay contains 5 rows (markers) and 4951 columns
#'   (cells); colData contains 4951 rows (cells) and 3 columns (metadata).
#' @seealso \code{\link{defined_image}} \code{\link{image_no_markers}}
"simulated_image"
