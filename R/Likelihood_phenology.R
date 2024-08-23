#' likelihood_phenology estimate likelihood for a set of parameters.
#' @title Estimate the likelihood of timeseries based on a set of parameters.
#' @author Marc Girondot
#' @return The likelihood of the data with the parameters
#' @param data Dataset generated with add_format
#' @param fixed.parameters Set of fixed parameters
#' @param fitted.parameters Set of parameters to be fitted
#' @param parallel If TRUE, parallel computing is used.
#' @param model_before The change of parameters before to estimate daily counts.
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
#' data_Gratiot <- add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg <- par_init(data_Gratiot, fixed.parameters=NULL)
#' # Estimate likelihood with this initial set of parameters
#' likelihood_phenology(data=data_Gratiot, fitted.parameters=parg, fixed.parameters=NULL)
#' # Or directly from a result object
#' likelihood_phenology(result=result_Gratiot)
#' # With new parametrization based on Omeyer et al. (2022)
#' # Omeyer, L. C. M., McKinley, T. J., Bréheret, N., Bal, G., Balchin, G. P., Bitsindou, A., 
#' # Chauvet, E., Collins, T., Curran, B. K., Formia, A., Girard, A., Girondot, M., Godley, B. J., 
#' # Mavoungou, J.-G., Poli, L., Tilley, D., VanLeeuwe, H. & Metcalfe, K. 2022. Missing data in 
#' # sea turtle population monitoring: a Bayesian statistical framework accounting for incomplete 
#' # sampling Front. Mar. Sci. (IF 3.661), 9, 817014.
#' 
#' parg <- c(tp=unname(parg["Peak"]), tf=unname(parg["Flat"]), 
#'           s1=unname(parg["LengthB"])/4.8, s2=unname(parg["LengthE"])/4.8, 
#'           alpha=unname(parg["Max_Complete"]), Theta=unname(parg["Theta"]))
#' likelihood_phenology(data=data_Gratiot, fitted.parameters=parg, fixed.parameters=NULL)
#' }
#' @export


likelihood_phenology <-  function(data=NULL                  , 
                                  fitted.parameters=NULL     , 
                                  fixed.parameters=NULL      , 
                                  parallel=TRUE              , 
                                  result=NULL                , 
                                  model_before=NULL          , 
                                  cofactors=NULL             , 
                                  add.cofactors=NULL         ,
                                  zero=1E-9                  , 
                                  out=TRUE                   ) {
    
    # data=NULL; fitted.parameters=NULL; fixed.parameters=NULL; parallel=TRUE; result=NULL; zero=1E-9; cofactors=NULL; add.cofactors=NULL; out=TRUE
    
    # if result est donné, on prend les données dedans et on remplace celles introduites en plus
    
    if (!is.null(result)) {
      if (!inherits(result, "phenology")) {
        stop("The object result must be the result of a fit_phenology()")
      }
      if (is.null(model_before)) {model_before <- result$model_before}
      if (is.null(data)) {data <- result$data}
      if (is.null(fitted.parameters)) {fitted.parameters <- result$par}
      if (is.null(fixed.parameters)) {fixed.parameters <- result$fixed.parameters}
      if (is.null(cofactors)) cofactors <- result$cofactors
      if (is.null(add.cofactors)) add.cofactors <- result$add.cofactors
    }
    
    
    if ((!is.null(add.cofactors)) & (!is.null(cofactors))) {
      cf1 <- cofactors
      for (ff in add.cofactors)
        cf1[, ff] <- cofactors[, ff] - mean(cofactors[, ff])
      cofactors <- cf1
    }
    
    
    # if (is.null(fixed.parameters)) {fixed.parameters <- NA}
    

    
    LnL <- getFromNamespace(".Lnegbin", ns="phenology")(x=fitted.parameters, 
                                                        pt=list(data=data, fixed=fixed.parameters, 
                                                                parallel=parallel, 
                                                                model_before=model_before, 
                                                                out=out, 
                                                                namespar=names(fitted.parameters), 
                                                                method_Snbinom=result$method_Snbinom, 
                                                                zero=zero, cofactors=cofactors, 
                                                                add.cofactors=add.cofactors))
    
    return(LnL)
    
  }
