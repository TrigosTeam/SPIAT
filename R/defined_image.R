#' SCE object of a simulated image with defined cell types based on marker
#' combinations.
#'
#' A dataset that contains a formatted sce object  \code{\link{simulated_image}}
#' with cell ids, cell positions, phenotypes, defined cell types in metadata and
#' marker intensities in assays.
#'
#' @format An sce object. Assay contains 5 rows (markers) and 4951 columns
#'   (cells); colData contains 4951 rows (cells) and 4 columns (metadata).
#' @seealso \code{\link{simulated_image}} \code{\link{image_no_markers}}
"defined_image"
