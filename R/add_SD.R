#' add_SD adds SD for a fixed parameter.
#' @title Add SD for a fixed parameter.
#' @author Marc Girondot
#' @return The parameters set with the new SD value
#' @param fixed.parameters Set of fixed parameters
#' @param parameters Set of current parameters
#' @param SD Standard deviation value to be added
#' @description This function is used to add standard deviation for a fixed parameter.
#' @examples
#' library(phenology)
#' # Generate a set of fixed parameter: Flat and Min
#' pfixed<-c(Flat=0, Min=0)
#'	# Add SD for the Flat parameter
#' pfixed<-add_SD(fixed.parameters=pfixed, parameters="Flat", SD=5)

#' @export


add_SD <-
function(fixed.parameters=NULL, parameters=NULL, SD=NULL) {

if (!is.null(fixed.parameters) && !is.null(SD) && !is.null(parameters)) {
	c<-SD
	names(c)<-paste("sd#",parameters,sep="")
	return(c(fixed.parameters, c))
}
}
