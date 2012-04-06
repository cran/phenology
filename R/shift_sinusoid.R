#' shift_sinusoid shift sinusoid information.
#' @title Shift sinusoid information.
#' @author Marc Girondot
#' @return Return a set of modified parameters
#' @param parameters set of parameters
#' @param from The number of series to change
#' @param to The number of series to change
#' @param help If TRUE, an help is displayed
#' @description This function is used to shift sinusoid parameters from '', '1' or '2'.
#' @export


shift_sinusoid <-
function(parameters=NULL, from="", to="1", help=FALSE) {
if(help) {
	cat("This function is used to roll sinusoid parameters from '', '1' or '2'.\n")
	cat("The syntax is par<-shift_sinusoid(parameters=value, from='1', to='2')\n")
	cat("or\n")
	cat("par<-shift_sinusoid(parameters=value, from='', to='1')\n")

} else {
if (!is.null(parameters) && (from!=to)) {

	level<-from
	varfrom<-c(paste("Phi", level, sep=""), paste("Delta", level, sep=""),paste("Alpha", level, sep=""),paste("Beta", level, sep=""),paste("Tau", level, sep=""))
	level<-to
	varto<-c(paste("Phi", level, sep=""), paste("Delta", level, sep=""),paste("Alpha", level, sep=""),paste("Beta", level, sep=""),paste("Tau", level, sep=""))
	
	for (i in 1:5) {
		if (length(names(parameters)[names(parameters)==varfrom[i]])==1) {
			names(parameters)[names(parameters)==varfrom[i]]<-varto[i]
		}
	}
	return(parameters)

}
}
}
