#' compute_gradient
#'
#' @description The function sweeps over circles of a range of radii surrounding
#'   reference cells and calculates the metrics at the radii. Metrics used with
#'   function need two conditions: 1) have a `radius` parameter. 2) return a
#'   single number. For metrics that do not return a single number, users can
#'   wrap them in a new function that returns a number and then pass the new
#'   function to `compute_gradient()`.
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param radii Numeric Vector specifying the range of radii for the metrics to
#'   be calculated.
#' @param FUN Variable name specifying the metric.
#' @param ... Arguments of FUN
#' @return A list of the metrics under all radii
#' @export
#'
#' @examples
#' gradient_positions <- c(30, 50, 100)
#' gradient_entropy <- compute_gradient(SPIAT::defined_image, 
#' radii = gradient_positions, FUN = calculate_entropy,  
#' cell_types_of_interest = c("Immune1","Immune2"), 
#' feature_colname = "Cell.Type")

compute_gradient <- function(spe_object, radii, FUN, ...){
    list.metric <- list() 
    for (i in seq_len(length(radii))){
        metric <- FUN(spe_object,radius = radii[i], ...)
        list.metric[[i]] <- metric 
    }
    return(list.metric)
}
