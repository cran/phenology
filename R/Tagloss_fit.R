#' Tagloss_fit fits a model of tag loss using a CMR database.
#' @title fit a model of tag loss using a CMR database.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
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
#' @param mc.cores Number of cores to use for parallel computing
#' @param groups Number of groups for parallel computing
#' @param hessian Does the hessian matrix should be estimated
#' @description This function fits a model of tag loss using a CMR database.\cr
#' The names of parameters can be:\cr
#' Model Pfaller et al. (2019):\cr
#' \describe{
#'   \item{Left tag lost when 2 are present}{\code{D1_L2}, \code{D2D1_L2}, \code{D3D2_L2}, \code{A_L2}, \code{B_L2}, \code{C_L2}, \code{delta_L2}}
#'   \item{Right tag lost when 2 are present}{\code{D1_R2}, \code{D2D1_R2}, \code{D3D2_R2}, \code{A_R2}, \code{B_R2}, \code{C_R2}, \code{delta_R2}}
#'   \item{Left tag lost when 1 is present}{\code{D1_L1}, \code{D2D1_L1}, \code{D3D2_L1}, \code{A_L1}, \code{B_L1}, \code{C_L1}, \code{delta_L1}}
#'   \item{Right tag lost when 1 is present}{\code{D1_R1}, \code{D2D1_R1}, \code{D3D2_R1}, \code{A_R1}, \code{B_R1}, \code{C_R1}, \code{delta_R1}}
#'   \item{One tag lost when 2 are present}{\code{D1_2}, \code{D2D1_2}, \code{D3D2_2}, \code{A_2}, \code{B_2}, \code{C_2}, \code{delta_2}}
#'   \item{One tag lost when 1 is present}{\code{D1_1}, \code{D2D1_1}, \code{D3D2_1}, \code{A_1}, \code{B_1}, \code{C_1}, \code{delta_1}}
#' }
#' \code{pA}, \code{pB} and \code{pC} are the daily probabilities of tag loss with 
#' \code{pA=-logit(A)}, \code{pB=-logit(B)} and \code{pC=-logit(C)} .\cr
#' \code{delta} is used as: \code{p = p + delta}. Nothe that \code{delta} can be negative\cr
#' Tag loss rate is \code{pA} at day 1\cr
#' Tag loss rate changes gradually from \code{pA} to \code{pB} that is reached at day \code{D1}\cr
#' Tag loss rate is \code{pB} from day \code{D1} to day \code{D2=D1+D2D1}\cr
#' Tag loss rate changes gradually from \code{pB} to \code{pC} that is reached at day \code{D3=D2+D3D2}\cr
#' 
#' When parameters from Rivalan et al. (2005) are used:\cr
#' \describe{
#' \item{One tag lost when 2 are present}{\code{a0_2}, \code{a1_2}, \code{a2_2}, \code{a3_2}, \code{a4_2}, \code{delta_2}}
#' \item{One tag lost when 1 is present}{\code{a0_1}, \code{a1_1}, \code{a2_1}, \code{a3_1}, \code{a4_1}, \code{delta_1}}
#' }
#' 
#' When parameters from Casale et al. (2017) are used:\cr
#' Model I\cr
#' \describe{
#' \item{One tag lost when 2 are present}{\code{CasaleModelIc_2}}
#' \item{One tag lost when 1 is present}{\code{CasaleModelIc_1}}
#' }
#' Model II\cr
#' \describe{
#' \item{One tag lost when 2 are present}{\code{CasaleModelIIa0_2}, \code{CasaleModelIIa1_2}, \code{CasaleModelIIa4_2}}
#' \item{One tag lost when 1 is present}{\code{CasaleModelIIa0_1}, \code{CasaleModelIIa1_1}, \code{CasaleModelIIa4_1}}
#' }
#' Model III\cr
#' \describe{
#' \item{One tag lost when 2 are present}{\code{CasaleModelIIIa0_2}, \code{CasaleModelIIIa1_2}, \code{CasaleModelIIIa4_2}}
#' \item{One tag lost when 1 is present}{\code{CasaleModelIIIa0_1}, \code{CasaleModelIIIa1_1}, \code{CasaleModelIIIa4_1}}
#' }
#' Model IV\cr
#' \describe{
#' \item{One tag lost when 2 are present}{\code{CasaleModelIVa0_2}, \code{CasaleModelIVa1_2}, \code{CasaleModelIVa2_2}, \code{CasaleModelIVa3_2}, \code{CasaleModelIVa4_2}}
#' \item{One tag lost when 1 is present}{\code{CasaleModelIVa0_1}, \code{CasaleModelIVa1_1}, \code{CasaleModelIVa2_1}, \code{CasaleModelIVa3_1}, \code{CasaleModelIVa4_1}}
#' }
#' Model V\cr
#'\describe{
#' \item{One tag lost when 2 are present}{\code{CasaleModelVa0_2}, \code{CasaleModelVa1_2}, \code{CasaleModelVa2_2}, \code{CasaleModelVa3_2}, \code{CasaleModelVa4_2}}
#' \item{One tag lost when 1 is present}{\code{CasaleModelVa0_1}, \code{CasaleModelVa1_1}, \code{CasaleModelVa2_1}, \code{CasaleModelVa3_1}, \code{CasaleModelVa4_1}}
#' }
#' 
#' If only one parameter is fitted, method must be "Brent" and \code{upper} and \code{lower} 
#' parameters must be set up with finite values.
#' 
#' model_before can be ""par['a0_1']=par['a0_2'];par['a1_1']=par['a1_2']".
#' model_after can be "p1=p2"
#' @family Model of Tag-loss
#' @references Rivalan, P., Godfrey, M.H., Prévot-Julliard, A.-C., Girondot, M., 2005. Maximum likelihood estimates of tag loss in leatherback sea turtles. Journal of Wildlife Management 69, 540-548.
#' @references Casale, P., Freggi, D., Salvemini, P., 2017. Tag loss is a minor limiting factor in sea turtle tagging programs relying on distant tag returns: the case of Mediterranean loggerhead sea turtles. European Journal of Wildlife Research 63.
#' @references Pfaller JB, Williams KL, Frick MG, Shamblin BM, Nairn CJ, Girondot M (2019) Genetic determination of tag loss dynamics in nesting loggerhead turtles: A new chapter in “the tag loss problem”. Marine Biology 166: 97 doi 10.1007/s00227-019-3545-x
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' data_f_21 <- Tagloss_format(outLR, model="21")
#' 
#' # model fitted by Rivalan et al. 2005
#' par <- c(a0_2=-5.43E-2, a1_2=-103.52, a4_2=5.62E-4, 
#'          delta_1=3.2E-4)
#' pfixed <- c(a2_2=0, a3_2=0, a2_1=0, a3_1=0)
#' model_before <- "par['a0_1']=par['a0_2'];par['a1_1']=par['a1_2'];par['a4_1']=par['a4_2']"
#' o <- Tagloss_fit(data=data_f_21, fitted.parameters=par, fixed.parameters=pfixed, 
#'                  model_before=model_before)
#' plot(o, t=1:1000, model="cumul")
#' plot(o, t=1:1000, model="1")
#' plot(o, t=1:1000, model="2", add=TRUE, col="red")
#' 
#' # Same data fitted with new model
#' par <- c(D1_1 = 100.15324837975547, A_1 = 5.9576927964120188, 
#'          B_1 = 8.769924225871069, B_2 = 8.2353860179664125)
#' pfixed <- c(D2D1_1 = 2568, D3D2_1 = 2568, D2D1_2 = 2568, D3D2_2 = 2568)
#' o_4p_p1p2 <- Tagloss_fit(data=data_f_21, fitted.parameters = par, 
#'                          fixed.parameters = pfixed, 
#'                          model_before = "par['C_1']=par['B_1'];
#'                          par['A_2']=par['A_1'];
#'                          par['C_2']=par['B_2'];
#'                          par['D1_2']=par['D1_1']", hessian=TRUE)
#'                          
#' # Without the N20 the computing is much faster
#' data_f_21_fast <- subset(data_f_21, subset=(is.na(data_f_21$N20)))
#' par <- c('D1_2' = 49.78891736351531, 
#'          'D2D1_2' = 1059.3635769732305, 
#'          'D3D2_2' = 12.434313273804602, 
#'          'A_2' = 5.2238379144659683, 
#'          'B_2' = 8.0050044071275543, 
#'          'C_2' = 8.4317863609499675, 
#'          'D1_1' = 701.80273287212935, 
#'          'D2D1_1' = 0.010951749100596819, 
#'          'D3D2_1' = 3773.6290607434876, 
#'          'A_1' = 205.42435592344776, 
#'          'B_1' = 9.9598342503239863, 
#'          'C_1' = 6.7234868237164722)
#' o <- Tagloss_fit(data=data_f_21_fast, fitted.parameters=par, hessian = TRUE)
#' plot(o, model="1", col="red")
#' plot(o, model="2", col="blue", add=TRUE)
#' legend("topright", legend=c("2->1", "1->0"), lty=1, col=c("blue", "red"))
#' }
#' @export

Tagloss_fit <- function (data = stop("A database formated using Tagloss_format() must be used"), 
                         fitted.parameters = NULL, fixed.parameters = NULL, model_before = NULL, 
                         model_after = NULL, control = list(trace = 1, maxit = 10000), 
                         method="Nelder-Mead", lower = -Inf, upper = Inf,
                         hessian = FALSE, mc.cores = detectCores(all.tests = FALSE, logical = TRUE), groups=NULL) {
  
  # data = NULL 
  # fitted.parameters = NULL; fixed.parameters = NULL; model_before = NULL 
  # model_after = NULL; control = list(trace = 1, maxit = 10000)
  # method="Nelder-Mead"; lower = -Inf; upper = Inf
  # hessian = FALSE; mc.cores = detectCores(all.tests = FALSE, logical = TRUE); groups=NULL
  
  if (!inherits(data, "TaglossData")) {
    stop("'data' must be a database formated using Tagloss_format().")
  }
  
  # Plus nécessaire
  if (FALSE) {
    data <- addS3Class(data, c("TaglossData", "data.frame"))
  }
  
  days.maximum <- Tagloss_daymax(data)
  
  # Tagloss_L(individuals = data, 
  #           par = fitted.parameters, 
  #           days.maximum = days.maximum, fixed.par = fixed.parameters, 
  #           model_before = model_before, 
  #           model_after = model_after)
  o <- optim(par = fitted.parameters, fn = Tagloss_L, individuals = data, 
             days.maximum = days.maximum, fixed.par = fixed.parameters, 
             model_before = model_before, 
             model_after = model_after, control = control, hessian = hessian, 
             method= method, lower = lower, upper = upper,
             mc.cores=mc.cores, groups=groups, names.par=names(fitted.parameters))
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
  o <- addS3Class(o, "Tagloss")
  
  return(o)
}
