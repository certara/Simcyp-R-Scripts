#' Create a vector of a rainbow of colors
#'
#' \code{rainbow} is useful for getting a set of colors for graphs when you need
#' a LOT of colors; there are prettier and more easily distinguishable options
#' if you only have a couple of comparisons.
#'
#' @param ncolors number of colors desired
#'
#' @return a character vector of colors
#' @export
#'
#' @examples
#' # Create a rainbow of 10 colors.
#' rainbow(10)
#' 
#' # Try using scales::show_col() to visualize the colors, ex:
#' scales::show_col(rainbow(10))
#' 
rainbow <- function(ncolors){
   
   if(ncolors <= 7){
      colRainbow <- colorRampPalette(c("#B91E02FF", "#F66D19FF", "#B8F735FF",
                                       "#22C4E3FF", "#434DB9FF"))
      
   } else {
      colRainbow <- colorRampPalette(c("gray20", "antiquewhite4", "firebrick3",
                                       "darkorange", "green3", "seagreen3",
                                       "cadetblue", "dodgerblue3", "royalblue4",
                                       "darkorchid4"))
      
   }
   
   return(colRainbow(ncolors))
}



#' Create a vector of blues fading into greens
#'
#' \code{blueGreens} is useful for getting a set of pretty blues that fade into
#' greens for graphs.
#'
#' @param ncolors number of colors desired
#'
#' @return a character vector of colors
#' @export
#'
#' @examples
#' # Create a set of 10 blues to greens.
#' blueGreens(10)
#'
#' # Try using scales::show_col() to visualize the colors, ex:
#' scales::show_col(blueGreens(10))
#' 
blueGreens <- function(ncolors){
   
   colblueGreens <- colorRampPalette(c("royalblue4", "dodgerblue3",
                                       "cadetblue", "seagreen3", "green3"))
   return(colblueGreens(ncolors))
}


#' Create a vector of blues
#'
#' \code{blues} is useful for getting a set of pretty blues for graphs. These
#' blues come from the 4th through 9th blues from grDevices::blues9, just to
#' give credit where it's due. 
#'
#' @param ncolors number of colors desired
#'
#' @return a character vector of colors
#' @export
#'
#' @examples
#' # Create a set of 10 blues.
#' blues(10)
#'
#' # Try using scales::show_col() to visualize the colors, ex:
#' scales::show_col(blues(10))
#' 
blues <- function(ncolors){
   
   colBlues <- colorRampPalette(c("#9ECAE1", "#6BAED6", "#4292C6", "#2171B5",
                                  "#08519C", "#08306B"))
   
   return(colBlues(ncolors))
   
}

