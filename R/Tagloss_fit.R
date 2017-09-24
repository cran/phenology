#' Tagloss_fit fits a model of tag loss using a CMR database.
#' @title fit a model of tag loss using a CMR database.
#' @author Marc Girondot
#' @return Return a list object with the model describing tag loss.
#' @param data An object formated using Tagloss_format
#' @param fitted.parameters Set of parameters to be fitted
#' @param fixed.parameters Set of fixed parameters
#' @param model_before Transformation of parameters before to use Tagloss_model()
#' @param model_after Transformation of parameters after to use Tagloss_model()
#' @param control Control parameters to be send to optim()
#' @param method optim() method
#' @param lower Lower value for parameter when Brent method is used
#' @param upper Upper value for parameter when Brent method is used
#' @param cores Number of cores to use for parallel computing
#' @param groups Number of groups for parallel computing
#' @param hessian Does the hessian matrix should be estimated
#' @description This function fits a model of tag loss using a CMR database.\cr
#' The names of parameters can be:\cr
#' Left tag lost when 2 are present: D1_L2, D2D1_L2, D3D2_L2, A_L2, B_L2, C_L2, delta_L2\cr
#' Right tag lost when 2 are present: D1_R2, D2D1_R2, D3D2_R2, A_R2, B_R2, C_R2, delta_R2\cr
#' Left tag lost when 1 is present: D1_L1, D2D1_L1, D3D2_L1, A_L1, B_L1, C_L1, delta_L1\cr
#' Right tag lost when 1 is present: D1_R1, D2D1_R1, D3D2_R1, A_R1, B_R1, C_R1, delta_R1\cr
#' One tag lost when 2 are present: D1_2, D2D1_2, D3D2_2, A_2, B_2, C_2, delta_2\cr
#' One tag lost when 1 is present: D1_1, D2D1_1, D3D2_1, A_1, B_1, C_1, delta_1\cr
#' A, B and C are -logit(pA), -logit(pB) and -logit(pC) of the 
#' corresponding daily probabilities (p) of tag loss.\cr
#' delta is used as: p = p * invlogit(-delta)\cr
#' The use of delta parameter is complicated.\cr
#' Tag loss rate is pA at day 1\cr
#' Tag loss rate changes gradually from pA to pB that is reached at day D1\cr
#' Tag loss rate is pB from day D1 to day D2=D1+D2D1\cr
#' Tag loss rate changes gradually from pB to pC that is reached at day D2+D3D2\cr
#' If only one parameter is fitted, method must be "Brent" and upper and lower 
#' parameters must be set up with finite values.
#' @family Model of Tag-loss
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' data_f_21 <- Tagloss_format(outLR, model="21")
#' # Without the N20 the computing is much faster
#' data_f_21_fast <- subset(data_f_21, subset=(is.na(data_f_21$N20)))
#' par <- structure(c(49.5658922243074, 808.136085362158, 106.283783786853, 
#' 5.22150592456511, 8.00608716525864, 8.32718202233396, 150.612916258503, 
#' 715.865805125223, 2242.06574225966, 119.212383120678, 10.1860735529433, 
#' 7.14231725937626), .Names = c("D1_2", "D2D1_2", "D3D2_2", "A_2", 
#' "B_2", "C_2", "D1_1", "D2D1_1", "D3D2_1", "A_1", "B_1", "C_1"))
#' o <- Tagloss_fit(data=data_f_21_fast, fitted.parameters=par)
#' plot(o, model="1", col="red")
#' plot(o, model="2", col="blue", add=TRUE)
#' legend("topright", legend=c("2->1", "1->0"), lty=1, col=c("blue", "red"))
#' }
#' @export

Tagloss_fit <- function (data = stop("A database formated using Tagloss_format() must be used"), 
          fitted.parameters = NULL, fixed.parameters = NULL, model_before = NULL, 
          model_after = NULL, control = list(trace = 1, maxit = 10000), 
          method="Nelder-Mead", lower = -Inf, upper = Inf,
          hessian = FALSE, cores = detectCores(all.tests = FALSE, logical = TRUE), groups=NULL) 
{
  days.maximum <- Tagloss_daymax(data)
  o <- optim(par = fitted.parameters, fn = Tagloss_L, individuals = data, 
             days.maximum = days.maximum, fixed.par = fixed.parameters, model_before = model_before, 
             model_after = model_after, control = control, hessian = hessian, 
             method= method, lower = lower, upper = upper,
             cores=cores, groups=groups, names.par=names(fitted.parameters))
  o$data <- data
  o$fixed.par <- fixed.parameters
  names(o$par) <- names(fitted.parameters)
  if (!is.null(o$hessian)) 
    o$SE <- SEfromHessian(o$hessian)
  o$AIC <- 2 * o$value + 2 * length(o$par)
  o$BIC <- 2 * o$value + log(nrow(data)) * length(o$par)
  o$AICc <- o$AIC + (2 * (length(o$par)) * (length(o$par) - 
                                              1))/(nrow(data) - length(o$par) - 1)
  o$model_before <- model_before
  o$model_after <- model_after
  o$days.maximum <- days.maximum
  class(o) <- "Tagloss"
  return(o)
}
