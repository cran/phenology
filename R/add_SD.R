#' add_SD adds SD for a fixed parameter.
#' @title Add SD for a fixed parameter.
#' @author Marc Girondot
#' @return The parameters set with the new SD value
#' @param parametersfixed Set of fixed parameters
#' @param parameter Set of current parameters
#' @param SD Standard deviation value to be added
#' @param help If TRUE, an help is displayed
#' @description This function is used to add standard deviation for a fixed parameter.
#' @export


add_SD <-
function(parametersfixed=NULL, parameter=NULL, SD=NULL, help=FALSE) {
if(help) {
	cat("This function is used to add standard deviation for a fixed parameter.\n")
	cat("The syntax is parfixed<-add_SD(parametersfixed=NULL, parameter=name, SD=value)\n")

} else {
if (!is.null(parametersfixed) && !is.null(SD) && !is.null(parameter)) {
	c<-SD
	names(c)<-paste("sd#",parameter,sep="")
	return(c(parametersfixed, c))
}
}
}
