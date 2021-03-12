#' Create point pattern object from single cell experiment object
#' 
#' Input - Single cell experiment object to be transformed into ppp object,
#'  marker phenotypes will be allocated as marks to ppp object
#' 
#' Output - ppp object, marked to reflect initial image
#' 
#' @param sce_object sce object to be transformed into ppp object


format_sce_to_ppp <- function(sce_object){
  # insert sce_object data into ppp object
  point_pattern <- ppp(colData(sce_object)[,2], colData(sce_object)[,3], c(min(colData(sce_object)[,2]), max(colData(sce_object)[,2])), c(min(colData(sce_object)[,3]), max(colData(sce_object)[,3])), marks = as.factor(colData(sce_object)[,1]))
  
  #return ppp object
  return(point_pattern)
}
