#' identify_bordering_cells
#'
#' @description Identify the cells bordering a group of cells of a particular
#'   phenotype.
#' @details The bordering cell detection algorithm is based on computing an
#'   alpha hull (Hemmer et al., 2020), a generalization of convex hull (Green
#'   and Silverman, 1979). The cells detected to be on the alpha hull are
#'   identified as the bordering cells.
#'
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param reference_cell String. Cells of this cell type will be used for border
#'   detection.
#' @param feature_colname String that speficies the column of `reference_cell`.
#' @param ahull_alpha Number specifying the parameter for the alpha hull
#'   algorithm. The larger the number, the more cells will be included in one
#'   cell cluster.
#' @param n_of_polygons Integer specifying the number of tumour regions defined
#'   by user.
#' @param draw Boolean if user chooses to manually draw the tumour area or not.
#'   Default is False.
#' @param n_to_exclude Integer. Clusters with cell count under this number will
#'   be deleted.
#' @export
#' @return A new SPE object is returned
#' @examples
#' spe_border <- identify_bordering_cells(SPIAT::defined_image, 
#' reference_cell = "Tumour", feature_colname = "Cell.Type", n_to_exclude = 10)

identify_bordering_cells <- function(spe_object, reference_cell, 
                                     feature_colname = "Cell.Type",
                                     ahull_alpha = NULL, n_of_polygons = 1, 
                                     draw = FALSE, n_to_exclude = 10){
    # CHECK
    if (is.null(SummarizedExperiment::colData(spe_object)[,feature_colname])){
        stop("Undefined column name. Please check if input contains `feature_colname`!")
    }
    if (!(reference_cell %in% 
          SummarizedExperiment::colData(spe_object)[,feature_colname])){
        stop("Reference cell not found!")
    }
    
    coords_df <- data.frame(SpatialExperiment::spatialCoords(spe_object))
    
    # plot #####
    phenotypes_of_interest <- c(reference_cell)
    colour_vector <- c("green")
    if (draw){
        graphics::par(xpd=TRUE)
        p <- plot_cell_basic(spe_object, phenotypes_of_interest, colour_vector,
                             feature_colname)
        graphics::par(xpd=FALSE)
    }
    
    ##### interactively draw boundaries ####
    if (n_of_polygons == 1 && !draw){
        l <- list()
        draw.polys <- data.frame("x" = c(max(coords_df$Cell.X.Position), 
                                         max(coords_df$Cell.X.Position),
                                         min(coords_df$Cell.X.Position), 
                                         min(coords_df$Cell.X.Position)),
                                 "y" = c(min(coords_df$Cell.Y.Position), 
                                         max(coords_df$Cell.Y.Position),
                                         max(coords_df$Cell.Y.Position), 
                                         min(coords_df$Cell.Y.Position)))
        poly <- sp::Polygon(draw.polys, hole = FALSE)
        l[[1]] <- poly
    }else{
        l <- list()
        for (i in seq_len(n_of_polygons)){
            draw.polys <- xROI::drawPolygon()
            poly <- sp::Polygon(draw.polys, hole = FALSE)
            l[[i]] <- poly}}
    
    polys <- sp::Polygons(l,ID = c("a"))
    sp_obj <- sp::SpatialPolygons(list(polys))
    
    # for loop, get the boundary cells and inside cells for each polygon #####
    data <- data.frame(SummarizedExperiment::colData(spe_object))
    data <- cbind(data,data.frame(SpatialExperiment::spatialCoords(spe_object)))
    data[,"Region"] <- "Outside"
    
    for (i in seq_len(n_of_polygons)){
        # get the coords
        buffered_polygon <- methods::slot(sp_obj@polygons[[1]]@Polygons[[i]],
                                          "coords")
        
        # identify the tumour cells in the drawn polygon
        inpolygon <- sp::point.in.polygon(
            data$Cell.X.Position, data$Cell.Y.Position, buffered_polygon[, 1], 
            buffered_polygon[, 2])
        allcells_in_polygon <- data[which(inpolygon!= 0),
                                   c("Phenotype","Cell.X.Position",
                                     "Cell.Y.Position",feature_colname)]
        tumour_in_polygon <- 
            allcells_in_polygon[which(allcells_in_polygon[,feature_colname] == 
                                          reference_cell),]
        tumour_in_polygon <- unique(tumour_in_polygon)
        
        # ahull of the tumour cells
        # define the value of alpha
        if (is.null(ahull_alpha)){
            n_cells <- dim(tumour_in_polygon)[1]
            if (n_cells<200){
                alpha <- 60} else if (n_cells > 5000){
                alpha <- 90} else {
                alpha <- (n_cells - 300)/160 + 60}
            methods::show(paste("The alpha of Polygon is:", alpha))
            ahull <- alphahull::ahull(tumour_in_polygon$Cell.X.Position,
                                     tumour_in_polygon$Cell.Y.Position, 
                                     alpha = alpha)
        }
        else {
            ahull <- alphahull::ahull(tumour_in_polygon$Cell.X.Position,
                                     tumour_in_polygon$Cell.Y.Position, 
                                     alpha = ahull_alpha)
        }
        # fix ahull
        ahull <- fix_ahull(ahull)
        # my polygon function
        xahull <- ahull$xahull
        arc <- ahull$arcs
        ahull_polygon <- get_polygon(xahull,arc,n_to_exclude)
        
        # identify the cells that compose the ahull
        border_ids <- c()
        for (i in c(seq_len(length(ahull_polygon)))){
            p <- data.frame(ahull_polygon[[i]])
            colnames(p) <- c("Cell.X.Position","Cell.Y.Position")
            
            # find the bordering cells in the original dataset to find the IDs
            common_cells <- dplyr::intersect(
                data[,c("Cell.X.Position","Cell.Y.Position")],
                p[,c("Cell.X.Position","Cell.Y.Position")])
            border_ids <- c(border_ids, rownames(common_cells))}
        
        # identify the cells that are in the ahull
        points_in_polygon <- data.frame()
        for (i in c(seq_len(length(ahull_polygon)))){
            p <- ahull_polygon[[i]]
            in_p <- sp::point.in.polygon(allcells_in_polygon$Cell.X.Position, 
                                         allcells_in_polygon$Cell.Y.Position, 
                                         p[,1], p[,2])
            points_in_polygon <- 
                unique(rbind(
                    points_in_polygon, 
                    allcells_in_polygon[which(in_p == 1),
                                        c("Phenotype", "Cell.X.Position",
                                          "Cell.Y.Position",
                                          feature_colname)]))}
        points_in_polygon_df <- as.data.frame(points_in_polygon)
        cells_in_boundary <- points_in_polygon_df
        
        in_border_ids <- rownames(cells_in_boundary)
        data[in_border_ids,"Region"] <- "Inside"
        data[border_ids,"Region"] <- "Border"
    }
    
    # plot and return #####
    SummarizedExperiment::colData(spe_object)$Region <- data[,"Region"]
    plot(data[which(data$Region=="Border"), 
              c("Cell.X.Position","Cell.Y.Position")], 
         pch = 19, cex = 0.3, main = paste(attr(spe_object, "name"),
                                           "tumour bordering cells"))
    return(spe_object)
}
