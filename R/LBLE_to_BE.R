#' LBLE_to_BE transforms a set of parameters from LengthB LengthE to Begin End format.
#' @title Transform a set of parameters from LengthB LengthE to Begin End format.
#' @author Marc Girondot
#' @return Return a set of modified parameters
#' @param parameters Set of current parameters
#' @param help If TRUE, an help is displayed
#' @description This function is used to transform a set of parameters
#' that uses LengthB, Peak and LengthE to a set of parameters
#' that uses Begin, Peak and End.

#' @export

LBLE_to_BE <-
function(parameters=NULL, help=FALSE) {
if(help) {
	cat("This function is used to transform a set of parameters\n")
	cat("that uses LengthB, Peak and LengthE to a set of parameters\n")
	cat("that uses Begin, Peak and End.\n")

} else {

if (!is.null(parameters)) {
	pk<-parameters["Peak"]
	lb<-parameters["LengthB"]
	le<-parameters["LengthE"]

	bg<-pk-lb
	en<-pk+le
	
	parameters[names(parameters)=="LengthB"]<-bg
	names(parameters)[names(parameters)=="LengthB"]<-"Begin"
	
	parameters[names(parameters)=="LengthE"]<-en
	names(parameters)[names(parameters)=="LengthE"]<-"End"
	
	return(parameters)
	
}
}
}
