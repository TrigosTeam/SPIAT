#' identify_bordering_cells
#'
#' @description Identify the cells bordering a group of cells of a particular phenotype
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param reference_cell Cells positive for this marker will be used as reference
#' @param draw Boolean if user chooses to draw the tumour area or not. Default is False.
#' @param n_of_polygons Number specifying the number of tumour regions defined by user
#' @param ahull_alpha Number specifying the ahull parameter. Larger number, more points included in the ahull.
#' @param column Column to select for phenotypes. Can be Phenotypes, Cell.Type, etc
#' @import SingleCellExperiment
#' @import alphahull
#' @importFrom xROI drawPolygon
#' @import sp
#' @importFrom dplyr intersect
#' @export

identify_bordering_cells <- function(sce_object, reference_cell, draw = F, 
                                     n_of_polygons = 1, ahull_alpha = NULL, 
                                     column = "Cell.Type"){
  # CHECK
  if (is.null(colData(sce_object)[,column])){
    stop("Please define the cell types!")
  }
  if (!(reference_cell %in% colData(sce_object)[,column])){
    stop("Reference cell not found!")
  }
  
  ##### plot #####
  phenotypes_of_interest <- c(reference_cell)
  colour_vector <- c("green")
  if (draw){
    par(xpd=TRUE)
    p <- plot_cell_basic(sce_object, phenotypes_of_interest, colour_vector, column = column)
    par(xpd=FALSE)
  }
  
  ##### interactively draw boundaries ####
  if (n_of_polygons == 1 && !draw){
    l <- list()
    draw.polys <- data.frame("x" = c(max(sce_object$Cell.X.Position), max(sce_object$Cell.X.Position),
                               min(sce_object$Cell.X.Position), min(sce_object$Cell.X.Position)),
                       "y" = c(min(sce_object$Cell.Y.Position), max(sce_object$Cell.Y.Position),
                               max(sce_object$Cell.Y.Position), min(sce_object$Cell.Y.Position)))
    poly <- Polygon(draw.polys, hole = FALSE)
    l[[1]] <- poly
  }
  else{
    l <- list()
    for (i in 1:n_of_polygons){
      draw.polys <- drawPolygon()
      poly <- Polygon(draw.polys, hole = FALSE)
      l[[i]] <- poly
    }
  }
  
  polys <- Polygons(l,ID = c("a"))
  sp <- SpatialPolygons(list(polys))
  
  ##### for loop, get the boundary cells and inside cells for each polygon #####
  data = data.frame(colData(sce_object))
  data[,"Region"] <- "Outside"
  
  for (i in 1:n_of_polygons){
    # get the coords
    buffered_polygon = slot(sp@polygons[[1]]@Polygons[[i]],"coords")
    
    # identify the tumour cells in the drawn polygon
    inpolygon = point.in.polygon(sce_object$Cell.X.Position, sce_object$Cell.Y.Position,
                                 buffered_polygon[, 1], buffered_polygon[, 2])
    allcells_in_polygon = data[which(inpolygon!= 0),c("Phenotype","Cell.X.Position",
                                                      "Cell.Y.Position",column)]
    tumour_in_polygon = allcells_in_polygon[which(allcells_in_polygon[,column] == reference_cell),]
    tumour_in_polygon = unique(tumour_in_polygon)
    
    # ahull of the tumour cells
    # define the value of alpha
    if (is.null(ahull_alpha)){
      n_cells = dim(tumour_in_polygon)[1]
      if (n_cells<200){
        alpha = 60
      } else if (n_cells > 5000){
        alpha = 90
      } else {
        alpha = (n_cells - 300)/160 + 60
      }
      print(paste("The alpha of Polygon is:", alpha))
      ahull = ahull(tumour_in_polygon$Cell.X.Position,
                    tumour_in_polygon$Cell.Y.Position, alpha = alpha)
    }
    else {
      ahull = ahull(tumour_in_polygon$Cell.X.Position,
                    tumour_in_polygon$Cell.Y.Position, alpha = ahull_alpha)
    }
    
    # identify the cells that compose the ahull
    cells_on_boundary = cbind(data.frame(ahull[["ashape.obj"]][["edges"]])$x1,data.frame(ahull[["ashape.obj"]][["edges"]])$y1)
    cells_on_boundary = data.frame(cells_on_boundary)
    colnames(cells_on_boundary) <- c("Cell.X.Position","Cell.Y.Position")
    
    # find the bordering cells in the original dataset to find the IDs
    common_cells <- dplyr::intersect(data[,c("Cell.X.Position","Cell.Y.Position")],
                                     cells_on_boundary[,c("Cell.X.Position","Cell.Y.Position")])
    
    border_ids <- rownames(common_cells)
    
    # fix ahull
    ahull <- fix_ahull(ahull)
    
    # my polygon function
    xahull <- ahull$xahull
    arc <- ahull$arcs
    ahull_polygon <- get_polygon(xahull,arc)
    points_in_polygon <- data.frame()
    
    # identify the cells that are in the ahull
    for (i in c(1:length(ahull_polygon))){
      p <- ahull_polygon[[i]]
      in_p <- point.in.polygon(allcells_in_polygon$Cell.X.Position, allcells_in_polygon$Cell.Y.Position, p[,1], p[,2])
      points_in_polygon <- unique(rbind(points_in_polygon,
                                        allcells_in_polygon[which(in_p == 1),c("Phenotype","Cell.X.Position","Cell.Y.Position",column)]))
      
    }
    
    tumour_in_polygon_df = as.data.frame(tumour_in_polygon)
    points_in_polygon_df = as.data.frame(points_in_polygon)
    cells_in_boundary = unique(rbind(tumour_in_polygon_df, points_in_polygon_df))
    
    in_border_ids <- rownames(cells_in_boundary)
    data[in_border_ids,"Region"] <- "Inside"
    data[border_ids,"Region"] <- "Border"
  }
  
  ##### plot and return #####
  colData(sce_object)$Region <- data[,"Region"]
  plot(data[which(data$Region=="Border"), c("Cell.X.Position","Cell.Y.Position")], 
       pch = 19, cex = 0.3, main = paste(attr(sce_object, "name"),"tumour bordering cells"))
  
  return(sce_object)
}