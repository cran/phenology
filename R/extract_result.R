#' extract_result get the fitted parameters from a result object.
#' @title Extract the set of parameters from a result object.
#' @author Marc Girondot
#' @return Return the set of fitted parameters
#' @param result A result file
#' @param help If TRUE, an help is displayed
#' @description The function "extract_result" permits to extract the set of parameters from a result object obtained after fit_phenology.
#' @export


extract_result <-
function(result=res, help=FALSE) {
if(help) {
	cat("This function is used to get the set of parameters\n")
	cat("from a result object obtained after fit_phenology.\n")

} else {
	return(result$par)
}
}
