#' logLik.phenology Return Log Likelihood of a fit
#' @title Return Log Likelihood of a fit generated by fit_phenology
#' @author Marc Girondot
#' @return The Log Likelihood value of the fitted model and data
#' @param object A result file generated by fit_phenology
#' @param ... Not used
#' @description Return Log Likelihood of a fit generated by fit_phenology
#' @family Phenology model
#' @examples
#' \dontrun{
#' library(phenology)
#' data(result_Gratiot)
#' logLik(result_Gratiot)
#' AIC(result_Gratiot)
#' }
#' @method logLik phenology
#' @export


logLik.phenology <- function(object, ...) {
  l <- -object$value
  n <- sum(unlist(lapply(object$data, function(x) sum(!is.na(x[,"nombre"])))))
  attributes(l) <- list(nall=n , nobs=n , df=length(object$par) , class="logLik")
  return(l)
}
