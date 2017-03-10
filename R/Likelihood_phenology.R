#' likelihood_phenology estimate likelihood for a set of parameters.
#' @title Estimate the likelihood of timeseries based on a set of parameters.
#' @author Marc Girondot
#' @return The likelihood of the data with the parameters
#' @param data Dataset generated with add_format
#' @param fixed.parameters Set of fixed parameters
#' @param fitted.parameters Set of parameters to be fitted
#' @param method_incertitude 2 [default] is the correct one from a statistical point of view; \cr
#'                           0 is an aproximate method more rapid; \cr
#'                           1 is an alternative more rapid but biased.
#' @param zero_counts example c(TRUE, TRUE, FALSE) indicates whether the zeros have 
#'                    been recorder for each of these timeseries. Defaut is TRUE for all.
#' @param infinite Number of iterations for dSnbinom() used for method_incertitude='sum'
#' @param cofactors data.frame with a column Date and a column for each cofactor
#' @param add.cofactors Names of the column of parameter cofactors to use as a cofactor
#' @param zero If the theoretical nest number is under this value, this value wll be used
#' @param result An object obtained after fit_phenology()
#' @description This function is used to estimate the likelihood based on a set of parameters.
#' @examples
#' \dontrun{
#' # Read a file with data
#' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
#' data(Gratiot)
#' # Generate a formated list nammed data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, fixed.parameters=NULL)
#' # Estimate likelihood with this initial set of parameters
#' likelihood_phenology(data=data_Gratiot, fitted.parameters=parg, fixed.parameters=NULL)
#' # Or directly from a result object
#' likelihood_phenology(result=result_Gratiot)
#' }
#' @export

likelihood_phenology <-
function(data=NULL, fitted.parameters=NULL, fixed.parameters=NULL, zero_counts=NULL, 
         method_incertitude=NULL, result=NULL, 
         cofactors=NULL, add.cofactors=NULL,
         infinite=200, zero=1E-9) {

 # data=NULL; fitted.parameters=NULL; fixed.parameters=NULL; zero_counts=NULL; method_incertitude=NULL; result=NULL; infinite=200; zero=1E-9
  
# if result est donné, on prend les données dedans et on remplace celles introduites en plus

if (!is.null(result)) {
  if (class(result) != "phenology") {
    stop("The object result must be the result of a fit_phenology()")
  }

  if (is.null(data)) {data <- result$data}
  if (is.null(fitted.parameters)) {fitted.parameters <- result$par}
  if (is.null(fixed.parameters)) {fixed.parameters <- result$fixed.parameters}
  if (is.null(zero_counts)) {zero_counts <- result$zero_counts}
  if (is.null(method_incertitude)) {method_incertitude <- result$method_incertitude}
  if (is.null(cofactors)) cofactors <- result$cofactors
  if (is.null(add.cofactors)) add.cofactors <- result$add.cofactors
}

# if (!is.null(data)) {
#   if (class(data) != "phenologydata") {
#     stop("The data object must be the result of a add_format()")
#   }
# }


# if (is.null(fixed.parameters)) {fixed.parameters <- NA}
if (is.null(method_incertitude)) {method_incertitude <- 0}
if (is.null(zero_counts)) {zero_counts <- TRUE}
	
	if (length(zero_counts)==1) {zero_counts <- rep(zero_counts, length(data))}
	if (length(zero_counts)!=length(data)) {
		stop("zero_counts parameter must be TRUE (the zeros are used for all timeseries) or FALSE (the zeros are not used for all timeseries) or possess the same number of logical values than the number of series analyzed.")
	}
	
LnL <- getFromNamespace(".Lnegbin", ns="phenology")(fitted.parameters, 
                                                    pt=list(data=data, fixed=fixed.parameters, 
                                                            incertitude=method_incertitude, zerocounts=zero_counts, 
                                                            infinite=infinite, out=TRUE, 
                                                            zero=zero, cofactors=cofactors, 
                                                            add.cofactors=add.cofactors))
	
	return(LnL)

}
