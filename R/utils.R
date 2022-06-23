# function to remove rows where intensity is NA
remove_intensity_na <- function(intensity_columns) {
    n_before <- nrow(intensity_columns)
    intensity_columns <- stats::na.omit(intensity_columns)
    n_after <- length(attributes(stats::na.omit(intensity_columns))$row.names)
    n <- n_before - n_after
    message(sprintf("Note: %i rows removed due to NA intensity values.",n))
    return(intensity_columns)
}

# convert spe object to a data frame with both colData and intensity matrix
bind_info <- function(spe_object){
    formatted_data <- data.frame(SummarizedExperiment::colData(spe_object))
    formatted_data <- cbind(formatted_data, 
                            data.frame(SpatialExperiment::spatialCoords(spe_object)))
    #convert rowname to column
    formatted_data <- formatted_data %>% tibble::rownames_to_column("Cell.ID") 
    
    # get the intensity matrix
    intensity_matrix <- SummarizedExperiment::assay(spe_object)
    markers <- rownames(intensity_matrix)
    cell_ids <- colnames(intensity_matrix)
    rownames(intensity_matrix) <- NULL
    colnames(intensity_matrix) <- NULL
    intensity_t <- data.frame(t(intensity_matrix))
    colnames(intensity_t) <- markers
    
    # bind
    formatted_data <- cbind(formatted_data, intensity_t)
    
    # delete column `sample_id`
    formatted_data$sample_id <- NULL
    
    return(formatted_data)
}

# convert spe object to a data frame with only colData
get_colData <- function(spe_object){
    formatted_data <- data.frame(SummarizedExperiment::colData(spe_object))
    formatted_data <- cbind(formatted_data, 
                            data.frame(SpatialExperiment::spatialCoords(spe_object)))
    formatted_data <- formatted_data %>% tibble::rownames_to_column("Cell.ID")
    
    # delete column `sample_id`
    formatted_data$sample_id <- NULL
    
    return(formatted_data)
}

# do not print out the messages
quiet_basic <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
} 

# fix the bug in `ahull` function from `alphahull` package
fix_ahull <- function(ahull){ # the order of cells returned by ahull is messy
    # this function reorders the cells
    arc <- ahull$arcs
    n_cells <- dim(arc)[1]
    # copy the arc ends
    ends <- arc
    ends <- rbind(ends, ends[1,])
    i <- 1
    while (i < n_cells){
        end2 <- ends[i,8]
        next_end1 <- ends[i+1,7]
        next_end2 <- ends[i+1,8]
        while (end2 != next_end1){ # the connection breaks here
            if (next_end2 == next_end1) {
                # if the next cell is a loner cell, or the current cell is the last cell,
                # the current cell does not have to be issue cell
                i <- i+1
                if (i >= (dim(arc)[1])) break
                next_end1 <- arc[i+1,7]
                next_end2 <- arc[i+1,8]
            }
            
            else{ # the next cell is not loner cell, confirm the current cell is the issue cell
                # break the current loop, go to next cell (skip loner cells)
                for (j in (i+1):n_cells){
                    t_end1 <- ends[j,7]
                    t_end2 <- ends[j,8]
                    # found the correct cell!
                    if (t_end1 == end2 || t_end2 == end2){ 
                        ends[n_cells+1,] <- ends[i+1,]
                        ends[i+1,] <- ends[j,]
                        ends[j,] <- ends[n_cells+1,]
                    }
                    if (t_end2 == end2){
                        ends[i+1,8] <- t_end1
                        ends[i+1,7] <- end2
                    }
                }
                break
            }
        }
        # next cell
        i <- i+1
    }
    ends <- ends[-(n_cells+1),]
    ahull$arcs <- ends
    return(ahull)
}

# gets the coordinates of points on ahull
get_polygon <- function(xahull, arc, n_to_exclude){ 
    df_list <- list()
    n <- 0
    s <- 1
    for (i in seq_len(dim(arc)[1]-1)){
        if (arc[i,8] != arc[i+1,7]){
            if (i-s > 4){
                n <- n+1
                df_list[[n]] <- arc[c(s:i),]
            }
            s <- i+1
        }
    }
    if (length(df_list) == 0) df_list[[1]] <- arc
    poly_list <- list()
    c <- 0 
    for (j in c(seq_len(length(df_list)))){
        df <- df_list[[j]]
        # eliminate the small tumour clusters (add later)
        if (dim(df)[1] <= n_to_exclude){
            next
        }
        else{
            c <- c+1
            cell_ID <- c()
            locs <- c()
            for (i in seq_len(dim(df)[1])){
                cell_ID <- df[i,7]
                locs <- rbind(locs,xahull[cell_ID,c(1,2)])
            }
            poly_list[[c]] <- locs
        }
    }
    return(poly_list)
}

# Produces a scatter plot of the cells in the tissue. Cells are coloured
# categorically by specified column. Cells not part of the celltypes of interest
# will be coloured "lightgrey"
plot_cell_basic <- function(spe_object, cell_types_of_interest, colour_vector, 
                             feature_colname, cex = 0.4, ...) {
    Cell.X.Position <- Cell.Y.Position <- NULL
    assign(feature_colname, NULL)
    formatted_data <- data.frame(SummarizedExperiment::colData(spe_object))
    formatted_data <- cbind(formatted_data, 
                            data.frame(SpatialExperiment::spatialCoords(spe_object)))
    # delete column `sample_id`
    formatted_data$sample_id <- NULL
    
    if (length(cell_types_of_interest) != length(colour_vector)) {
        stop("The colour vector is not the same length as the celltypes of interest")
    }
    
    real_celltypes <- cell_types_of_interest
    for (phenotype in cell_types_of_interest) {
        if (!(phenotype %in% unique(formatted_data[[feature_colname]]))) {
            methods::show(paste(phenotype, "cells were not found"), sep = "")
            real_celltypes <- real_celltypes[real_celltypes != phenotype]
        }
    }
    
    colour_vector <- colour_vector[match(real_celltypes, cell_types_of_interest)]
    
    cell_types_of_interest <- real_celltypes
    
    if (any(!formatted_data[[feature_colname]] %in% cell_types_of_interest)) {
        formatted_data[!formatted_data[[feature_colname]] %in% cell_types_of_interest,
        ][[feature_colname]] <- "OTHER"
    }
    
    formatted_data$color <- ""
    
    for (phenotype in cell_types_of_interest) {
        idx <- which(cell_types_of_interest == phenotype)
        formatted_data[formatted_data[[feature_colname]] == phenotype,
        ]$color <- colour_vector[idx]
    }
    
    if (any(formatted_data[[feature_colname]] == "OTHER")) {
        formatted_data[formatted_data[[feature_colname]] == "OTHER", ]$color <- "grey"
        all_phenotypes <- unique(c(cell_types_of_interest, "OTHER"))
        all_colours <- c(colour_vector, "grey")[seq_len(length(all_phenotypes))]
    }
    
    else {
        all_phenotypes <- cell_types_of_interest
        all_colours <- colour_vector
    }
    
    name_of_object <- attr(spe_object, "name")
    if (!is.null(name_of_object)) name_of_object <- paste("", name_of_object)
    
    plot(formatted_data$Cell.X.Position,formatted_data$Cell.Y.Position,
         pch = 19,cex = cex, col = formatted_data$color,
         xlab = "X Position", ylab = "Y Position",
         main = paste0("Plot", name_of_object, " by " ,feature_colname), ...)
    
    x.max <- max(formatted_data$Cell.X.Position)
    y.max <- max(formatted_data$Cell.Y.Position)
    x <- x.max + 2
    y <- y.max/2
    
    graphics::par(mar=c(5.1, 4.1, 4.1, 5.1), xpd=TRUE)
    graphics::legend(x,y, legend=all_phenotypes,
           col=all_colours, cex=0.6, pch = 19, box.lty=0, bg='gray')
    
}

# define a function to get the number of certain name under a certain column
count_category <- function(spe_object, cat, feature_colname){
    data <- data.frame(SummarizedExperiment::colData(spe_object))
    count_table <- table(data[feature_colname])
    count <- unname(count_table[match(cat, names(count_table))])
    return(count)
}

split_image <- function(spe_object, number_of_splits){
    Cell.X.Position <- Cell.Y.Position <- NULL
    #Reads the image file
    cell_loc <- get_colData(spe_object)
    #turn the spe object name into string
    image_name <- deparse(substitute(spe_object))
    #CHECK
    if(nrow(cell_loc) == 0) stop("There are no cells in the data.")
    #Selects x and y region to plot
    minX <- min(cell_loc$Cell.X.Position, na.rm = TRUE)
    maxX <- max(cell_loc$Cell.X.Position, na.rm = TRUE)
    minY <- min(cell_loc$Cell.Y.Position, na.rm = TRUE)
    maxY <- max(cell_loc$Cell.Y.Position, na.rm = TRUE)
    #Splits the range of x and y coordinates into n + 1 evenly spaced out lengths
    x_split <- seq(minX, maxX, length.out = number_of_splits + 1)
    y_split <- seq(minY, maxY, length.out = number_of_splits + 1)
    
    divided_image_obj <- list()
    for(y in seq_len(number_of_splits)){
        local_coor_y <- y_split[c(y+1, y)]
        #Round coordinates of cuts to nearest whole number for labeling
        rounded_coor_y <- round(local_coor_y)
        temp_cell_loc <- cell_loc
        
        for(x in seq_len(number_of_splits)){
            local_coor_x <- x_split[c(x+1, x)]
            #Rounds coordinates of cuts to nearest whole number for labeling
            rounded_coor_x <- round(local_coor_x)
            temp_cell_loc <- cell_loc

            #Locates all points in the data file fitting the given min/max parameters to be plotted i.e. splitting up the image.
            divided_image <- temp_cell_loc[min(local_coor_x) < temp_cell_loc$Cell.X.Position & 
                                               temp_cell_loc$Cell.X.Position <= max(local_coor_x)
                                           & min(local_coor_y) < temp_cell_loc$Cell.Y.Position & 
                                               temp_cell_loc$Cell.Y.Position <= max(local_coor_y), ]

            divided_image_obj[[paste(image_name,"r", x,"c", y, sep="")]] <- 
                divided_image
        }
    }
    return(divided_image_obj)
}
