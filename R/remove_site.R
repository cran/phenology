#' remove_site removes beach information from a set of parameters.
#' @title Removes site information from a set of parameters.
#' @author Marc Girondot
#' @return Return a set of modified parameters
#' @param parameters Set of parameters
#' @param help If TRUE, an help is displayed
#' @description This function is used to remove the information of the site
#' from a set of parameters. It can be used to other timeseries after.
#' @export



remove_site <-
function(parameters=NULL, help=FALSE) {
if(help) {
	cat("This function is used to remove the information of the site\n")
	cat("from a set of parameters. It can be used to other timeseries after.\n")

} else {
if (!is.null(parameters)) {
	for(i in 1:length(parameters)) {
		names(parameters)[i]<-strsplit(names(parameters[i]), "_")[[1]][1]
	}
	return(parameters)
}
}
}
