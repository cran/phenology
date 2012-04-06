#' Likelihood_phenology estimate likelihood for a set of parameters.
#' @title Estimate the likelihood of timeseries based on a set of parameters.
#' @author Marc Girondot
#' @return Return a list of formated data
#' @param data =dataset generated with add_format
#' @param parametersfixed Set of fixed parameters
#' @param parametersfit Set of parameters to be fitted
#' @param method_incertitude 2 [default] is the correct one from a statistical point of view; 
#'                           0 is an aproximate method more rapid; 
#'                           1 is an alternative more rapid but biased.
#' @param help If TRUE, an help is displayed
#' @description This function is used to estimate the likelihood based on a set of parameters.
#' @export


Likelihood_phenology <-
function(data=dta, parametersfit=x, parametersfixed=NULL, method_incertitude=0, help=FALSE) {

if(help) {
	cat("This function is used to estimate the likelihood based on\n")
	cat("a set of parameters.\n")
	cat("The general syntax is:\n")
	cat("Likelihood_phenology(data=dta, parametersfit=x, parametersfixed=parfixed, method_incertitude=0)\n")
} else {

if (is.null(parametersfixed)) {parametersfixed<-NA}

	.phenology.env <<- new.env()
	.phenology.env$data<<-data
	.phenology.env$fixed<<-parametersfixed
	.phenology.env$incertitude<<-method_incertitude


	LnL<-.Lnegbin(parametersfit)
	
	cat(paste("-LnL=", LnL, sep=""))

}
}
