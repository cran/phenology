#' likelihood_phenology estimate likelihood for a set of parameters.
#' @title Estimate the likelihood of timeseries based on a set of parameters.
#' @author Marc Girondot
#' @return The likelihood of the data with the parameters
#' @param data Dataset generated with add_format
#' @param fixed.parameters Set of fixed parameters
#' @param fitted.parameters Set of parameters to be fitted
#' @param zero_counts example c(TRUE, TRUE, FALSE) indicates whether the zeros have 
#'                    been recorder for each of these timeseries. Defaut is TRUE for all.
#' @param tol Tolerance of recurrence for dSnbinom() used for convolution of negative binomial distribution
#' @param parallel If TRUE, parallel computing is used.
#' @param cofactors data.frame with a column Date and a column for each cofactor
#' @param add.cofactors Names of the column of parameter cofactors to use as a cofactor
#' @param zero If the theoretical nest number is under this value, this value wll be used
#' @param result An object obtained after fit_phenology()
#' @param out If TRUE, return the global likelihood; if FALSE, the likelihood for each series
#' @description This function is used to estimate the likelihood based on a set of parameters.
#' @family Phenology model
#' @examples
#' \dontrun{
#' # Read a file with data
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
  function(data=NULL, 
           fitted.parameters=NULL, fixed.parameters=NULL, zero_counts=NULL, 
           parallel=TRUE, 
           result=NULL, 
           cofactors=NULL, add.cofactors=NULL,
           tol=1E-6, zero=1E-9, out=TRUE) {
    
    # data=NULL; fitted.parameters=NULL; fixed.parameters=NULL; zero_counts=NULL; parallel=TRUE; result=NULL; zero=1E-9; cofactors=NULL; add.cofactors=NULL; tol=1E-6; out=TRUE
    
    # if result est donné, on prend les données dedans et on remplace celles introduites en plus
    
    if (!is.null(result)) {
      if (class(result) != "phenology") {
        stop("The object result must be the result of a fit_phenology()")
      }
      
      if (is.null(data)) {data <- result$data}
      if (is.null(fitted.parameters)) {fitted.parameters <- result$par}
      if (is.null(fixed.parameters)) {fixed.parameters <- result$fixed.parameters}
      if (is.null(zero_counts)) {zero_counts <- result$zero_counts}
      if (is.null(cofactors)) cofactors <- result$cofactors
      if (is.null(add.cofactors)) add.cofactors <- result$add.cofactors
    }
    
    # if (!is.null(data)) {
    #   if (class(data) != "phenologydata") {
    #     stop("The data object must be the result of a add_format()")
    #   }
    # }
    
    if ((!is.null(add.cofactors)) & (!is.null(cofactors))) {
      cf1 <- cofactors
      for (ff in add.cofactors)
        cf1[, ff] <- cofactors[, ff] - mean(cofactors[, ff])
      cofactors <- cf1
    }
    
    
    # if (is.null(fixed.parameters)) {fixed.parameters <- NA}
    if (is.null(zero_counts)) {zero_counts <- TRUE}
    
    if (length(zero_counts)==1) {zero_counts <- rep(zero_counts, length(data))}
    if (length(zero_counts)!=length(data)) {
      stop("zero_counts parameter must be TRUE (the zeros are used for all timeseries) or FALSE (the zeros are not used for all timeseries) or possess the same number of logical values than the number of series analyzed.")
    }
    
    LnL <- getFromNamespace(".Lnegbin", ns="phenology")(x=fitted.parameters, 
                                                        pt=list(data=data, fixed=fixed.parameters, 
                                                                parallel=parallel, 
                                                                zerocounts=zero_counts, 
                                                                tol=tol, out=out, 
                                                                namespar=names(fitted.parameters), 
                                                                zero=zero, cofactors=cofactors, 
                                                                add.cofactors=add.cofactors))
    
    return(LnL)
    
  }
