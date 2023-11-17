
#' Calculate the Ki of an  inhibitor with the Cheng-Prusoff equation
#'
#' \code{calcKi} uses IC50, substrate concentration, and Km to calculate a Ki
#' value for an inhibitor
#'
#' @param IC50 The IC50 value for the inhibitor (uM)
#' @param S The concentration of the substrate used (uM)
#' @param Km The Michaelis constant for the affinity between the substrate and
#'   the enzyme (uM)
#'
#' @return Returns the inhibitory constant Ki (uM) (numeric)
#' @export
#'
#' @examples
#' calcKi(20, 5, 5)
#' # 10
#'

calcKi <- function(IC50, S, Km = NA){
    Ki <- IC50/(1+S/Km)
    return(Ki)
}
