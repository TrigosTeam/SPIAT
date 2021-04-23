#' identify_bordering_cells_interactive
#'
#' @description Identify the cells bordering a group of cells of a particular phenotype
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param reference_cell Cells positive for this marker will be used as reference
#' @param n_of_polygons Number specifying the number of tumour regions defined by user
#' @param buffer_width Number specifying the allowed margin of error
#' @param ahull_alpha Number specifying the ahull parameter. Larger number, more points included in the ahull.
#' @param large Boolean specifying if the image requires splitting
#' @import SingleCellExperiment
#' @import alphahull
#' @import xROI
#' @import raster
#' @import sp
#' @import dplyr
#' @import hull2spatial
#' @export

#sce_object = formatted_defined
#reference_cell = "MEL"
#n_of_polygons = 2
#buffer_width = 30
#ahull_alpha = 60

# colData() is in package 'SummarizedExperiment' but imported by SingleCellExperiment

identify_bordering_cells_interactive <- function(sce_object, reference_cell, n_of_polygons = 1, buffer_width = 1, ahull_alpha = NULL, large = FALSE){
  # CHECK
  if (is.null(sce_object$Cell.Type)){
    stop("Please define the cell types!")
  }
  if (!(reference_cell %in% sce_object$Cell.Type)){
    stop("Reference cell not found!")
  }
  
  ##### plot #####
  phenotypes_of_interest <- c(reference_cell)
  colour_vector <- c("green")
  par(xpd=TRUE)
  p <- plot_cell_basic(sce_object, phenotypes_of_interest, colour_vector, column = "Cell.Type")
  par(xpd=FALSE)
  
  ##### interactively draw boundaries ####
  if (n_of_polygons == 1){
    list <- list()
    draw <- data.frame("x" = c(max(sce_object$Cell.X.Position), max(sce_object$Cell.X.Position),
                               min(sce_object$Cell.X.Position), min(sce_object$Cell.X.Position)),
                       "y" = c(min(sce_object$Cell.Y.Position), max(sce_object$Cell.Y.Position),
                               max(sce_object$Cell.Y.Position), min(sce_object$Cell.Y.Position)))
    poly <- Polygon(draw, hole = FALSE)
    list[[1]] <- poly
  }
  
  else{
    list = list()
    for (i in 1:n_of_polygons){
      draw <- drawPolygon()
      poly <- Polygon(draw, hole = FALSE)
      list[[i]] <- poly
    }
  }
  polys <- Polygons(list,ID = c("a"))
  sp <- SpatialPolygons(list(polys))
  
  ##### for loop, get the boundary cells and inside cells for each polygon #####
  data = data.frame(colData(sce_object))
  data[,"Region"] <- "Out"
  buffer = buffer(sp, width = buffer_width)

  
  for (i in 1:n_of_polygons){
    # get the coords
    #print(buffer@polygons[[1]])
    #print(buffer@polygons[[1]]@Polygons[[i]])
    buffered_polygon = slot(buffer@polygons[[1]]@Polygons[[i]],"coords")
    # identify the tumour cells in the drawn polygon
    inpolygon = point.in.polygon(sce_object$Cell.X.Position, sce_object$Cell.Y.Position, 
                                 buffered_polygon[,"x"], buffered_polygon[,"y"])
    allcells_in_polygon = data[which(inpolygon!= 0),c("Phenotype","Cell.X.Position",
                                                      "Cell.Y.Position","Cell.Type")]
    tumour_in_polygon = allcells_in_polygon[which(allcells_in_polygon$Cell.Type == reference_cell),]
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
      ahull= ahull(tumour_in_polygon$Cell.X.Position, 
                   tumour_in_polygon$Cell.Y.Position, alpha = alpha)
    }
    else {
      ahull= ahull(tumour_in_polygon$Cell.X.Position, 
                   tumour_in_polygon$Cell.Y.Position, alpha = ahull_alpha)
    }
    
    
    # identify the cells that compose the ahull
    cells_on_boundary = cbind(data.frame(ahull$ashape.obj$edges)$x1,data.frame(ahull$ashape.obj$edges)$y1)
    cells_on_boundary = data.frame(cells_on_boundary)
    colnames(cells_on_boundary) <- c("Cell.X.Position","Cell.Y.Position")
    common_cells <- dplyr::intersect(data[,c("Cell.X.Position","Cell.Y.Position")], 
                                     cells_on_boundary[,c("Cell.X.Position","Cell.Y.Position")])

    border_ids <- rownames(common_cells)

    # change ahull into spatial lines
    ahull_line <- ahull2lines(ahull)
    ahull_polygon <- ahull_line@lines[[1]]@Lines[[1]]@coords
    
    # identify the cells that are in the ahull
    intumour = point.in.polygon(allcells_in_polygon$Cell.X.Position, allcells_in_polygon$Cell.Y.Position, ahull_polygon[,"x"], ahull_polygon[,"y"])
    points_in_polygon = allcells_in_polygon[which(intumour!= 0),c("Phenotype","Cell.X.Position","Cell.Y.Position","Cell.Type")]
    tumour_in_polygon_df = as.data.frame(tumour_in_polygon)
    points_in_polygon_df = as.data.frame(points_in_polygon)
    cells_in_boundary = rbind(tumour_in_polygon_df, points_in_polygon_df)
    cells_in_boundary = unique(cells_in_boundary)
    in_border_ids <- rownames(cells_in_boundary)
    data[in_border_ids,"Region"] <- "In"
    data[border_ids,"Region"] <- "Border"
  }
  
  
  ##### plot and return #####
  colData(sce_object)$Region <- data[,"Region"]
  plot(data[which(data$Region=="Border"), c("Cell.X.Position","Cell.Y.Position")], pch = 19, cex = 0.3)
  
  return(sce_object)
}

