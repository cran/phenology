#' Tagloss_mcmc Bayesian model of tag loss using a CMR database.
#' @title Bayesian model of tag loss using a CMR database.
#' @author Marc Girondot
#' @return Return a list object with the Bayesian model describing tag loss.
#' @param data An object formated using Tagloss_format
#' @param parameters A data.frame with priors; see description and examples
#' @param fixed.parameters Set of fixed parameters
#' @param model_before Transformation of parameters before to use Tagloss_model()
#' @param model_after Transformation of parameters after to use Tagloss_model()
#' @param cores Number of cores to use for parallel computing
#' @param groups Number of groups for parallel computing
#' @param n.iter Number of iterations for each chain
#' @param n.chains Number of chains
#' @param n.adapt Number of iteration to stabilize likelihood
#' @param thin Interval for thinning likelihoods
#' @param trace Or FALSE or period to show progress
#' @param traceML TRUE or FALSE to show ML
#' @param intermediate Or NULL of period to save intermediate result
#' @param filename Name of file in which intermediate results are saved
#' @param adaptive Should an adaptive process for SDProp be used
#' @param adaptive.lag  Lag to analyze the SDProp value in an adaptive context
#' @param adaptive.fun Function used to change the SDProp
#' @param previous The content of the file in which intermediate results are saved
#' @description This function fits a model of tag loss using a CMR database using Bayesian mcmc.\cr
#' The parameters must be stored in a data.frame with named rows for each parameter with the following columns:\cr
#' \itemize{
#'   \item Density. The density function name, example \code{dnorm}, \code{dlnorm}, \code{dunif}
#'   \item Prior1. The first parameter to send to the \code{Density} function
#'   \item Prior2. The second parameter to send to the \code{Density} function
#'   \item SDProp. The standard error from new proposition value of this parameter
#'   \item Min. The minimum value for this parameter
#'   \item Max. The maximum value for this parameter
#'   \item Init. The initial value for this parameter
#' }
#' @family Model of Tag-loss
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
#' pMCMC <- data.frame()
#' o <- Tagloss_mcmc(data=data_f_21, parameters=pMCMC, fixed.parameters=pfixed, 
#'                  model_before=model_before)
#' }
#' @export

Tagloss_mcmc <- function (data = stop("A database formated using Tagloss_format() must be used"), 
                          parameters = stop("Priors must be supplied"), 
                          fixed.parameters = NULL, 
                          model_before = NULL, 
                          model_after = NULL, 
                          cores = detectCores(all.tests = FALSE, logical = TRUE), 
                          groups = detectCores(all.tests = FALSE, logical = TRUE), 
                          n.iter=10000, n.chains = 1, n.adapt = 100, thin=30, 
                          trace=FALSE, traceML=FALSE, 
                          adaptive = FALSE, adaptive.lag = 500, 
                          adaptive.fun = function(x) {ifelse(x>0.234, 1.3, 0.7)},
                          intermediate=NULL, filename="intermediate.Rdata",
                          previous=NULL) 
{
  
  Tagloss_MCMC <- function(x, individuals=NULL, 
                           days.maximum = NULL, 
                           fixed.parameters = NULL, 
                           model_before = NULL, 
                           model_after = NULL, 
                           names.par = NULL, 
                           groups = groups, 
                           cores = cores) {
    
    Tagloss_L(individuals=individuals, par=x, 
              days.maximum = days.maximum, 
              fixed.parameters =fixed.parameters, 
              model_before = model_before, 
              model_after = model_after, 
              names.par = names.par, 
              groups = groups, 
              cores = cores, 
              progressbar = FALSE)
    
  }
  
  
  MCMC <- MHalgoGen(likelihood = Tagloss_MCMC, 
                    parameters = parameters, 
                    individuals=data, 
                    days.maximum = Tagloss_daymax(data), 
                    fixed.parameters = fixed.parameters, 
                    model_before = model_before, 
                    model_after = model_after, 
                    groups = groups, 
                    cores = cores, 
                    n.iter=n.iter, 
                    n.chains = n.chains, n.adapt = n.adapt, thin = thin, 
                    adaptive.lag = adaptive.lag, 
                    adaptive=adaptive, 
                    adaptive.fun = adaptive.fun, 
                    trace = trace, 
                    traceML=traceML, 
                    previous=previous)
  
  return(MCMC)
}
