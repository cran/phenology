#' logLik.Tagloss returns Log Likelihood of a fit for tag loss
#' @title Return Log Likelihood of a fit generated by Tagloss_fit
#' @author Marc Girondot
#' @return The Log Likelihood value for the fitted model with data
#' @param object A result file generated by Tagloss_fit
#' @param ... Not used
#' @description Return Log Likelihood of a fit generated by Tagloss_fit
#' @family Model of Tag-loss
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' data_f_21 <- Tagloss_format(outLR, model="21")
#' # Without the N20 the computing is much faster
#' data_f_21_fast <- subset(data_f_21, subset=(is.na(data_f_21$N20)))
#' par <- structure(c(48.8292784204825, 1039.02842229274, -89.3162940697861, 
#' 5.21817463244988, 8.00575451188548, 8.32971268127933, 161.265553603601, 
#' 602.935748681661, 2643.57415102633, 16.752815732218, 10.181616195839, 
#' 7.14279063312016), .Names = c("D1_2", "D2D1_2", "D3D2_2", "A_2", 
#' "B_2", "C_2", "D1_1", "D2D1_1", "D3D2_1", "A_1", "B_1", "C_1"))
#' o <- Tagloss_fit(data=data_f_21_fast, fitted.parameters=par)
#' logLik(o)
#' AIC(o)
#' }
#' @method logLik Tagloss
#' @export


logLik.Tagloss <- function(object, ...) {
  l <- -object$value
  n <- nrow(object$data)
  attributes(l) <- list(nall=n , nobs=n , df=length(object$par) , class="logLik")
  return(l)
}
