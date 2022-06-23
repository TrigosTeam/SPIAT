#' Average nearest neighbor index for point pattern (clustering or dispersion)
#'
#' @description Calculate the the average nearest neighbor (ANN) index of a
#'   specified type of cells. The index indicates the clustering effect of a
#'   point pattern. The pattern can be clustering, random or dispersion.
#' @details ANN index is a statistical test to test for the presence of clusters
#'   of cells, (Clark and Evans, 1954). The ANN index evaluates the spatial
#'   aggregation or dispersion effect of objects based on the average distances
#'   between pairs of the nearest objects and can be used to test for the
#'   clustering of specific cell types (e.g. immune or tumor cells). Next, the z
#'   score and p-value of the ANN index is calculated to validate the
#'   significance of the pattern.
#'
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param reference_celltypes String Vector. Cells with these cell types will be
#'   used for ANNI calculation.
#' @param feature_colname String. Specify the selected column for
#'   `reference_celltypes`.
#' @param p_val Numeric. The p value threshold to determine the significance of
#'   a pattern.
#' @export
#' @return A list with the pattern type and a p value
#' @examples
#' average_nearest_neighbor_index(SPIAT::defined_image, reference_celltypes = 
#' "Tumour", feature_colname = "Cell.Type")

average_nearest_neighbor_index <- function(spe_object, reference_celltypes, 
                                           feature_colname, p_val = 5e-6){
    
    ppp <- format_spe_to_ppp(spe_object)
    formatted_data <- get_colData(spe_object)
    
    data <- formatted_data[,c(feature_colname,"Cell.X.Position",
                              "Cell.Y.Position")]
    data <- data[which(data[,feature_colname] %in% reference_celltypes),
                 c("Cell.X.Position","Cell.Y.Position") ]
    
    if(nrow(data) == 0){
        methods::show("No reference cells found")
        ANN_index <- list(pattern=NA,`p-value`=NA)
    }else{
        object<-format_colData_to_spe(data)
        
        ppp_object <- format_spe_to_ppp(object)
        ann.p <- mean(spatstat.geom::nndist(ppp_object, k=1))
        
        n <- ppp_object$n # Number of points
        
        x <- ppp$window$xrange[2] - ppp$window$xrange[1]
        y <- ppp$window$yrange[2] - ppp$window$yrange[1]
        area <- x*y
        ann.e <- 0.5/sqrt(n/area)
        se <- 0.26136/sqrt(n*n/area)
        z <- (ann.p - ann.e)/se
        p <- stats::pnorm(-abs(z))
        
        if (p <= p_val){
            if (ann.p <= ann.e){
                pattern <- "Clustered"
            }
            else{pattern <- "Dispersed"}
        }
        else{pattern <- "Random"}
        
        ANN_index <- list(pattern, p)
        names(ANN_index) <- c("pattern" ,"p-value")
    }
    return(ANN_index)
}
