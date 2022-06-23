#' Format SPE object as a ppp object (`spatstat` package)
#'
#' @description Formats an spe object into a ppp object which has the x,y
#'   coordinates, phenotypes as markers window specifies the range of x and y
#'   coordinates
#'
#' @export
#' @param spe_object SpatialExperiment object in the form of the output of
#'   format_image_to_spe.
#' @param window_pol Optional Boolean Specifying if the window is polygon.
#' @param feature_colname String specifying the feature column of interest.
#' @return A ppp object is returned (defined in `spatstat` package)
#' @examples
#' ppp_object<-format_spe_to_ppp(SPIAT::defined_image, 
#' feature_colname = "Cell.Type")

format_spe_to_ppp <- function(spe_object, window_pol = FALSE, 
                              feature_colname="Phenotype") {
    
    # get x, y coordinates and phenotypes from spe object
    spe_data <- get_colData(spe_object)
    spe_data <- spe_data[!duplicated(spe_data[,c("Cell.X.Position", 
                                                 "Cell.Y.Position")]),]
    x <- spe_data$Cell.X.Position
    y <- spe_data$Cell.Y.Position
    marks <- spe_data[[feature_colname]]
    
    # get windows
    x_summary <- summary(spe_data$Cell.X.Position)
    x_min <- as.numeric(x_summary[1])
    x_max <- as.numeric(x_summary[6])
    y_summary <- summary(spe_data$Cell.Y.Position)
    y_min <- as.numeric(y_summary[1])
    y_max <- as.numeric(y_summary[6])
    
    if (window_pol == TRUE){
        
        # get ploy window
        X <- data.frame(x,y)
        hpts <- grDevices::chull(X)
        poly_window <- list(x=rev(X[hpts, 1]), y=rev(X[hpts, 2]))
        
        # format spe to ppp
        ppp_object <- spatstat.geom::ppp(x, y, poly = poly_window, 
                                         marks = marks)
    }
    
    else{  
        ppp_object <- spatstat.geom::ppp(
            x, y, window = spatstat.geom::owin(c(x_min, x_max),c(y_min, y_max)),
            marks = marks)}
    return(ppp_object)
}
