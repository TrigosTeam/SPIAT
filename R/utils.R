# function to remove rows where intensity is NA
remove_intensity_na <- function(intensity_columns) {
    n_before <- nrow(intensity_columns)
    intensity_columns <- stats::na.omit(intensity_columns)
    n_after <- length(attributes(stats::na.omit(intensity_columns))$row.names)
    n <- n_before - n_after
    message(sprintf("Note: %i rows removed due to NA intensity values.",n))
    return(intensity_columns)
}

# convert sce object to a dataframe with both colData and intensity matrix
bind_colData_intensity <- function(sce_object){
    formatted_data <- data.frame(SummarizedExperiment::colData(sce_object))
    #convert rowname to column
    formatted_data <- formatted_data %>% tibble::rownames_to_column("Cell.ID") 
    
    # get the intensity matrix
    intensity_matrix <- SummarizedExperiment::assay(sce_object)
    markers <- rownames(intensity_matrix)
    cell_ids <- colnames(intensity_matrix)
    rownames(intensity_matrix) <- NULL
    colnames(intensity_matrix) <- NULL
    intensity_t <- data.frame(t(intensity_matrix))
    colnames(intensity_t) <- markers
    
    # bind
    formatted_data <- cbind(formatted_data, intensity_t)
    
    return(formatted_data)
}

# convert sce object to a dataframe with only colData
get_colData <- function(sce_object){
    formatted_data <- data.frame(SummarizedExperiment::colData(sce_object))
    formatted_data <- formatted_data %>% tibble::rownames_to_column("Cell.ID")
    return(formatted_data)
}

# do not print out the messages
quiet_basic <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
} 

# fix the bug in ahull function from alphahull package
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
            cell_ID = c()
            locs <- c()
            for (i in seq_len(dim(df)[1])){
                cell_ID = df[i,7]
                locs <- rbind(locs,xahull[cell_ID,c(1,2)])
            }
            poly_list[[c]] <- locs
        }
    }
    return(poly_list)
}

# Produces a scatter plot of the cells in the tissue. Cells are coloured
# categorically by specified column. Cells not part of the celltypes of interest will be coloured "lightgrey"
plot_cell_basic <- function (sce_object, cell_types_of_interest, colour_vector, 
                             feature_colname, cex = 0.4) {
    Cell.X.Position <- Cell.Y.Position <- NULL
    assign(feature_colname, NULL)
    formatted_data <- data.frame(colData(sce_object))
    
    if (length(cell_types_of_interest) != length(colour_vector)) {
        stop("The colour vector is not the same length as the celltypes of interest")
    }
    
    real_celltypes <- cell_types_of_interest
    for (phenotype in cell_types_of_interest) {
        if (!(phenotype %in% unique(formatted_data[[feature_colname]]))) {
            print(paste(phenotype, "cells were not found"), sep = "")
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
        all_phenotypes <- c(cell_types_of_interest, "OTHER")
        all_colours <- c(colour_vector, "grey")
    }
    
    else {
        all_phenotypes <- cell_types_of_interest
        all_colours <- colour_vector
    }
    
    name_of_object <- attr(sce_object, "name")
    
    plot(formatted_data$Cell.X.Position,formatted_data$Cell.Y.Position,
         pch = 19,cex = cex, col = formatted_data$color,
         xlab = "X Position", ylab = "Y Position",
         main = paste("Plot", name_of_object, "by" ,feature_colname))
    
    x.max <- max(sce_object$Cell.X.Position)
    y.max <- max(sce_object$Cell.Y.Position)
    x <- x.max
    y <- y.max/2
    graphics::par(mar=c(5.1, 4.1, 4.1, 5.1), xpd=TRUE)
    graphics::legend(x,y, legend=all_phenotypes,
           col=all_colours, cex=0.6, pch = 19,box.lty=0, bg='gray')
    
}

# define a function to get the number of certain name under a certain column
count_category <- function(sce_object, cat, feature_colname){
    data <- data.frame(SummarizedExperiment::colData(sce_object))
    count_table <- table(data[feature_colname])
    count <- unname(count_table[match(cat, names(count_table))])
    return(count)
}
