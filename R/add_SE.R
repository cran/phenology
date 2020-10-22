#' add_SE adds standard error for a fixed parameter.
#' @title Add standard error for a fixed parameter.
#' @author Marc Girondot
#' @return The parameters set with the new SE value
#' @param fixed.parameters Set of fixed parameters
#' @param parameters Set of current parameters
#' @param SE Standard error value to be added
#' @description This function is used to add standard error for a fixed parameter.
#' @family Phenology model
#' @examples
#' library(phenology)
#' # Generate a set of fixed parameter: Flat and Min
#' pfixed<-c(Flat=0, Min=0)
#'	# Add SE for the Flat parameter
#' pfixed<-add_SE(fixed.parameters=pfixed, parameters="Flat", SE=5)

#' @export


add_SE <-
function(fixed.parameters=NULL, parameters=NULL, SE=NULL) {

if (!is.null(fixed.parameters) && !is.null(SE) && !is.null(parameters)) {
	c<-SE
	names(c)<-paste("se#",parameters,sep="")
	return(c(fixed.parameters, c))
}
}
