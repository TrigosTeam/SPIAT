#' SPE object of a simulated image with defined cell types based on marker
#' combinations.
#'
#' A dataset that contains a formatted spe object with cell ids, phenotypes,
#' defined cell types in `colData()` and marker intensities in `assays()`. (The
#' cell locations are the same with the cells in \code{\link{simulated_image}}).
#'
#' @format An spe object. Assay contains 5 rows (markers) and 4951 columns
#'   (cells); colData contains 4951 rows (cells) and 3 columns (features).
#' @seealso \code{\link{simulated_image}} \code{\link{image_no_markers}}
"defined_image"
