#' colhex2col
#'
#' @description Converts a set of hexadecimal colours into normal human readable colours
#'
#' @param colhex - A vector of hexadecimal colours.
#' @return the human readable colour converted last from \code{colhex}.
#' @importFrom grDevices col2rgb colors

#HELPER FUNCTION
#Title: colhex2col.R
#Author: sklarz-bgu
#Date: 2nd December, 2019
#Availability: https://gist.github.com/sklarz-bgu/01a550f59cdf5bc85a48e15f5e94a6ba
colhex2col <- function(colhex) {
    
    # setting this to NULL as otherwise get "no visible binding for global variable" in R check
    . <- NULL
    
    # Convert hex to RGB
    mycol   <- colhex %>% col2rgb()
    # Convert all x11 colors to RGB, adn transform
    colors()                 %>%        # Get X11 colors (hex)
        col2rgb              %>%        # Convert to RGB matrix
        data.frame           %>%        # Convert to data.frame
        setNames(.,colors()) %>%        # Set color names
        t                    %>%        # Transform so colors are in rows (columns: R,G,B)
        data.frame           %>%        # Re-convert to data.frame
        apply(.,1,function(x) sum(abs(x-mycol)) ) %>%  # For each color, calc the sum of diff between mycol RGB and the color RGB
        sort                 %>%        # Sort so color with smallest diff comes first
        '['(1)               %>%        # Get the first color in the df = closest
        names                %>%        # Return the name of the color
        return
}
