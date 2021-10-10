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
        if (i >= (dim(arc)[1])){
          break
        }
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

get_polygon <- function(xahull, arc){ # this function gets the coordinates that on ahull
  df_list <- list()
  n <- 0
  s <- 1
  for (i in 1: (dim(arc)[1]-1)){
    if (arc[i,8] != arc[i+1,7]){
      if (i-s > 5){
        n <- n+1
        df_list[[n]] <- arc[c(s:i),]
      }
      s <- i+1
    }
  }
  if (length(df_list) == 0) df_list[[1]] <- arc
  poly_list <- list()
  c <- 0 
  for (j in c(1:length(df_list))){
    df <- df_list[[j]]
    c <- c+1
    cell_ID = c()
    locs <- c()
    for (i in 1:(dim(df)[1])){
      cell_ID = df[i,7]
      locs <- rbind(locs,xahull[cell_ID,c(1,2)])
    }
    poly_list[[c]] <- locs
    polygon(locs)
    
  }
  return(poly_list)
}
