#' phenology runs a shiny application for basic functions of phenology package
#' @title Run a shiny application for basic functions of phenology package
#' @author Marc Girondot
#' @return Nothing
#' @description Run a shiny application for basic functions of phenology package.
#' Thanks to Adriana Cortés Gomés and Joana Hancock for their help with translation.
#' @examples
#' \dontrun{
#' library(phenology)
#' phenology()
#' }
#' @export


phenology <- function() {
  
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package is absent; Please install it first")
  }
  
getFromNamespace("runApp", ns="shiny")(appDir = system.file("shiny", package="phenology"), 
                                       launch.browser =TRUE)

}
