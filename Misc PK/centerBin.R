#' Cut a set of numeric data into bins based on the middle value of the bin
#'
#' This function cuts a vector of numeric data such that the breaks are set
#' using the middle value of a given bin rather than setting upper or lower
#' boundaries for each bin. For example, say you've got a vector of blood draw
#' times that weren't drawn at \emph{exactly} the nominal time but were close to
#' it. Set the bins here to the nominal times, and the data will be filed into
#' the slot that's closest to the nominal time at which that sample should have
#' been drawn. For a similar option, see \code{\link{cutNumeric}} or see the
#' base R function \code{\link[base]{cut}}.
#'
#' @param x A numeric string
#' @param breaks the desired breaks in that string set by what the middle value
#'   should be. Any value of \code{x} exactly in the middle of two breaks will
#'   be assigned to the larger break, i.e., the function rounds up.
#' @examples
#' x <- c(48, 39, 28, 24, 5, 64, 133, 51, 59, 92, NA, 39)
#' centerBin(x, breaks = seq(0, 340, 5))
#' centerBin(x, breaks = seq(0, 340, 10))
#'
#' # Let's look at an example where the times in a PK study weren't perfect.
#' ActualTime <- c(-5, 17, 32, 65, 118)
#' NominalTime <- c(0, 15, 30, 60, 120)
#' centerBin(ActualTime, breaks = NominalTime)
#'
#'
#' @return Returns a vector of numeric data
#' @export
#'


centerBin <- function(x, breaks){

      breaks <- sort(unique(breaks))

      CenterBins <- data.frame(LfBin = c(NA, breaks),
                               RtBin = c(breaks, NA),
                               stringsAsFactors = FALSE)
      CenterBins$LfBin[1] <- breaks[1]
      CenterBins$RtBin[length(breaks) + 1] <- breaks[length(breaks)]
      CenterBins$Lower <- (CenterBins$RtBin - CenterBins$LfBin)/2 +
            CenterBins$LfBin
      CenterBins$Upper <- c(CenterBins$Lower[-1], NA)
      CenterBins$Upper[nrow(CenterBins)] <-
            max(c(max(x, na.rm = TRUE), max(breaks, na.rm = TRUE)))
      CenterBins <- dplyr::arrange(CenterBins, Lower)

      Outbin <- x
      for(i in seq_along(x)){
            Assignment <- CenterBins$RtBin[which(x[i] >= CenterBins$Lower &
                                                       x[i] <= CenterBins$Upper)]
            Assignment <- ifelse(length(Assignment) > 1,
                                 Assignment[2], Assignment[1])

            if(complete.cases(x[i])){
                  if(x[i] < min(breaks, na.rm = TRUE)){
                        Assignment <- min(breaks, na.rm = TRUE)
                  }

                  if(x[i] > max(breaks, na.rm = TRUE)){
                        Assignment <- max(breaks, na.rm = TRUE)
                  }

            } else {
                  Assignment <- NA
            }

            Outbin[i] <- Assignment

            rm(Assignment)
      }

      return(Outbin)
}
