#' Import spatstat library and create multitype point pattern based on 
#' single cell experiment object
#' 
#' input - single cell experiment object from SPIAT
#' output - multitype point pattern representative of sce_object
#' 

import_spatstat <- function(sce_object) {
  #import spatstat package
  library(spatstat)
  
  #create point pattern object
  point_pattern <- ppp(colData(sce_object)[,2], colData(sce_object)[,3], c(min(colData(sce_object)[,2]), max(colData(sce_object)[,2])), c(min(colData(sce_object)[,3]), max(colData(sce_object)[,3])), marks = as.factor(colData(sce_object)[,1]))
  return(point_pattern)
}

