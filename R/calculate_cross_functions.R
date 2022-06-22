#' calculate_cross_functions
#'
#' @description Compute and plot the cross functions between two specified cell
#'   types. This function implements the cross functions from [spatstat] package.
#'
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param method String that is the method for dependence calculation. Options:
#'   "Gcross", "Kcross", "Kcross.inhom", "Lcross", "Jcross". Default method is
#'   "Kcross".
#' @param cell_types_of_interest String Vector. Cell types of interest.
#' @param feature_colname String that is the name of the column of the types.
#' @param plot_results Boolean. TRUE if result to be plotted, FALSE if not. In
#'   either case, an object with the results is returned
#' @param dist Number (OPTIONAL) The largest distance between two cell types at
#'   which K function is evaluated. If NULL, use the default distances set by
#'   cross functions.
#' @importFrom spatstat.core Gcross Kcross.inhom Lcross Jcross Kcross
#' @importFrom spatstat.geom ppp
#' @return An object of class "fv" defined in `spatstat` package.
#' @export
#' @examples
#' df_cross <- calculate_cross_functions(SPIAT::defined_image, 
#' method = "Kcross", cell_types_of_interest = c("Tumour","Immune3"),
#' feature_colname ="Cell.Type", dist = 100)

calculate_cross_functions <- function(spe_object, method = "Kcross", 
                                      cell_types_of_interest, feature_colname, 
                                      plot_results = TRUE, dist = NULL) {
    #CHECK
    formatted_data <-get_colData(spe_object)
    if (!all(cell_types_of_interest %in% formatted_data[[feature_colname]])) {
        stop("Cell type not found!")
    }
    
    # format spe to ppp object
    ppp_object <- format_spe_to_ppp(spe_object, 
                                    feature_colname = feature_colname)
    ppp_object$marks <- as.factor(ppp_object$marks)
    
    # r
    if (is.null(dist)) r <- NULL
    else r <- seq(0, dist, length.out = 100)
    
    if (method == "Gcross"){
        p <- Gcross(ppp_object, cell_types_of_interest[1],
                    cell_types_of_interest[2],correction = "border", r = r)
        if(plot_results){
            plot(p, main = paste("cross G function",attr(spe_object,"name")))
        }
    }
    else if (method == "Kcross"){
        p <- Kcross(ppp_object, cell_types_of_interest[1],
                    cell_types_of_interest[2],correction = "border", r = r)
        if(plot_results){
            if (is.null(dist)) plot(p, main = paste("cross K function",
                                                    attr(spe_object,"name")))
            else plot(p, main = paste("cross K function",
                                      attr(spe_object,"name")), 
                      xlim = c(0,dist))
        }
    }
    else if (method == "Kcross.inhom"){
        p <- Kcross.inhom(ppp_object, cell_types_of_interest[1],
                          cell_types_of_interest[2],correction = "border", 
                          r = r)
        if(plot_results){
            if (is.null(dist)) plot(p, main = paste("cross K function",
                                                    attr(spe_object,"name")))
            else plot(p, main = paste("cross K function",
                                      attr(spe_object,"name")), 
                      xlim = c(0,dist))
        }
    }
    else if (method == "Lcross"){
        p <- Lcross(ppp_object, cell_types_of_interest[1],
                    cell_types_of_interest[2],correction = "border", r = r)
        if(plot_results){
            plot(p, main = paste("cross L function",attr(spe_object,"name")))
        }
    }
    else if (method == "Jcross"){
        p <- Jcross(ppp_object, cell_types_of_interest[1],
                    cell_types_of_interest[2],correction = "border", r = r)
        if(plot_results){
            plot(p, main = paste("cross J function",attr(spe_object,"name")))
        }
    }
    return(p)
}
