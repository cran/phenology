#' likelihood_phenology estimate likelihood for a set of parameters.
#' @title Estimate the likelihood of timeseries based on a set of parameters.
#' @author Marc Girondot
#' @return The likelihood of the data with the parameters
#' @param data Dataset generated with add_format
#' @param parametersfixed Set of fixed parameters
#' @param parametersfit Set of parameters to be fitted
#' @param method_incertitude 2 [default] is the correct one from a statistical point of view; \cr
#'                           0 is an aproximate method more rapid; \cr
#'                           1 is an alternative more rapid but biased.
#' @param zero_counts example c(TRUE, TRUE, FALSE) indicates whether the zeros have 
#'                    been recorder for each of these timeseries. Defaut is TRUE for all.
#' @param help If TRUE, an help is displayed
#' @description This function is used to estimate the likelihood based on a set of parameters.
#' @examples 
#' # Read a file with data
#' # Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", , header=FALSE)
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot<-add_format(origin=NULL, add=Gratiot, name="Complete", reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Estimate likelihood with this initial set of parameters
#' likelihood_phenology(data=data_Gratiot, parametersfit=parg, parametersfixed=NULL)
#' @export




likelihood_phenology <-
function(data=NULL, parametersfit=NULL, parametersfixed=NULL, zero_counts=TRUE, method_incertitude=0, help=FALSE) {

if(help) {
	cat("This function is used to estimate the likelihood based on\n")
	cat("a set of parameters.\n")
	cat("The general syntax is:\n")
	cat("likelihood_phenology(data=dta, parametersfit=x, parametersfixed=parfixed, \n")
	cat("+         zero_counts=TRUE, method_incertitude=0)\n")
} else {

if (is.null(parametersfixed)) {parametersfixed<-NA}

.phenology.env<- NULL
rm(.phenology.env)

	.phenology.env <<- new.env(parent=.GlobalEnv)
	assign("data", data, envir = as.environment(.phenology.env))
	assign("fixed", parametersfixed, envir = as.environment(.phenology.env))
	assign("incertitude", method_incertitude, envir = as.environment(.phenology.env))
	
	if (length(zero_counts)==1) {zero_counts<-rep(zero_counts, length(data))}
	if (length(zero_counts)!=length(data)) {
		print("zero_counts parameter must be TRUE (the zeros are used for all timeseries) or FALSE (the zeros are not used for all timeseries) or possess the same number of logical values than the number of series analyzed.")
		return()
	}
	
	assign("zerocounts", zero_counts, envir = as.environment(.phenology.env))

	LnL<-.Lnegbin(parametersfit)
	
	return(LnL)

}
}
