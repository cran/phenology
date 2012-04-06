#' BE_to_LBLE transforms a set of parameters from Begin End format to LengthB LengthE.
#' @title Transform a set of parameters from Begin End format to LengthB LengthE.
#' @author Marc Girondot
#' @return Return the set of modified parameters
#' @param parameters Set of current parameters
#' @param help If TRUE, an help is displayed
#' @description This function is used to transform a set of parameters 
#' that uses Begin, Peak and End to a set of parameters 
#' that uses LengthB, Peak and LengthE.
#' @export

BE_to_LBLE <-
function(parameters=NULL, help=FALSE) {
if(help) {
	cat("This function is used to transform a set of parameters\n")
	cat("that uses Begin, Peak and End to a set of parameters\n")
	cat("that uses LengthB, Peak and LengthE.\n")

} else {

if (!is.null(parameters)) {
	pk<-parameters["Peak"]
	bg<-parameters["Begin"]
	en<-parameters["End"]

	lb<-pk-bg
	le<-en-pk
	
	parameters[names(parameters)=="Begin"]<-lb
	names(parameters)[names(parameters)=="Begin"]<-"LengthB"
	
	parameters[names(parameters)=="End"]<-le
	names(parameters)[names(parameters)=="End"]<-"LengthE"
	
	return(parameters)
	
}
}
}
