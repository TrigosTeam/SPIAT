#' gg_color_hue
#'
#' @description Generates a vector of colours in hexadecimal
#'
#' @param number_of_colours - The number of hexadecimal colours to generate.
#' @return A vector of hexadecimal colours
#' @import RColorBrewer
#' @import grDevices

#HELPER FUNCTION
#Title: gg_color_hue.R
#Author: tinyheero
#Date: 2nd December, 2019
#Availability: https://github.com/tinyheero/tinyutils/blob/master/R/gg_color_hue.R
gg_color_hue <- function(number_of_colours) {
  hues = seq(15, 375, length = number_of_colours + 1)
  hcl(h = hues, l = 65, c = 100)[1:number_of_colours]
}
