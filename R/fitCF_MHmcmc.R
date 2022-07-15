#' fitCF_MHmcmc runs the Metropolis-Hastings algorithm for ECFOCF (Bayesian MCMC)
#' @title Run the Metropolis-Hastings algorithm for ECFOCF data
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A list with resultMCMC being mcmc.list object, resultLnL being likelihoods and parametersMCMC being the parameters used
#' @param n.iter Number of iterations for each step
#' @param parametersMCMC A set of parameters used as initial point for searching with information on priors
#' @param result An object obtained after a SearchR fit
#' @param n.chains Number of replicates
#' @param n.adapt Number of iterations before to store outputs
#' @param thin Number of iterations between each stored output
#' @param adaptive Should an adaptive process for SDProp be used
#' @param adaptive.lag  Lag to analyze the SDProp value in an adaptive content
#' @param adaptive.fun Function used to change the SDProp
#' @param trace TRUE or FALSE or period, shows progress
#' @param traceML TRUE or FALSE to show ML
#' @param filename If intermediate is not NULL, save intermediate result in this file
#' @param intermediate Period for saving intermediate result, NULL for no save
#' @param previous Previous result to be continued. Can be the filename in which intermediate results are saved.
#' @family Model of Clutch Frequency
#' @description Run the Metropolis-Hastings algorithm for RMU.data.\cr
#' The number of iterations is n.iter+n.adapt+1 because the initial likelihood is also displayed.\cr
#' I recommend thin=1 because the method to estimate SE uses resampling.\cr
#' As initial point is maximum likelihood, n.adapt = 0 is a good solution.\cr
#' The parameters intermediate and filename are used to save intermediate results every 'intermediate' iterations (for example 1000). Results are saved in a file of name filename.\cr
#' The parameter previous is used to indicate the list that has been save using the parameters intermediate and filename. It permits to continue a mcmc search.\cr
#' These options are used to prevent the consequences of computer crash or if the run is very very long and computer processes at time limited.\cr
#' @examples 
#' \dontrun{
#' library("phenology")
#' data(MarineTurtles_2002)
#' ECFOCF_2002 <- TableECFOCF(MarineTurtles_2002)
#' 
#' # Paraetric model for clutch frequency
#' o_mu1p1_CFp <- fitCF(x = c(mu = 2.1653229641404539, 
#'                  sd = 1.1465246643327098, 
#'                  p = 0.25785366120357966), 
#'                  fixed.parameters=NULL, 
#'                  data=ECFOCF_2002, hessian = TRUE)
#'                            
#' pMCMC <- fitCF_MHmcmc_p(result=o_mu1p1_CFp, accept=TRUE)
#' fitCF_MCMC <- fitCF_MHmcmc(result = o_mu1p1_CFp, n.iter = 1000, 
#'                            parametersMCMC = pMCMC, n.chains = 1, n.adapt = 0, 
#'                            adaptive=TRUE, 
#'                            thin = 1, trace = TRUE)
#'                            
#' plot(fitCF_MCMC, parameters="mu")
#' plot(fitCF_MCMC, parameters="sd")
#' plot(fitCF_MCMC, parameters="p", xlim=c(0, 0.5), breaks=seq(from=0, to=0.5, by=0.05))
#' plot(fitCF_MCMC, parameters="p", transform = invlogit, xlim=c(0, 1), 
#'      breaks=c(seq(from=0, to=1, by=0.05)))
#' 
#' }
#' @export


fitCF_MHmcmc <- function(result=stop("An output from fitCF() must be provided"), 
                          n.iter=10000, 
                          parametersMCMC=stop("A parameter set from fitCF_MHmcmc_p() must be provided"), 
                          n.chains = 1,
                          n.adapt = 0, thin=1, 
                          adaptive=FALSE, 
                          adaptive.lag=500, 
                          adaptive.fun=function(x) {ifelse(x>0.234, 1.3, 0.7)},
                          trace=FALSE, 
                          traceML=FALSE, 
                          intermediate=NULL, filename="intermediate.Rdata", previous=NULL) {
  
  if (is.character(previous)) {
    itr <- NULL
    load(previous)
    previous <- itr
    rm(itr)
    print("Continue previous mcmc run")
  }

  if (!inherits(result, "ECFOCF")) {
    stop("An output from fitCF() must be provided")
  }

  fun <- lnLCF

  print(parametersMCMC)
  
  # x, fixed, model.trend, RMU.data, colname.year=NULL, RMU.names=NULL, index=NULL
  # pt <- list(fixed=result$fixed.parameters, RMU.data=result$RMU.data, model.trend=result$model.trend, colname.year=result$colname.year, RMU.names=result$RMU.names)

out <- MHalgoGen(n.iter=n.iter, parameters=parametersMCMC, n.chains = n.chains, n.adapt = n.adapt, 
                 thin=thin, trace=trace, traceML=traceML, likelihood=fun,
                 adaptive = adaptive, adaptive.fun = adaptive.fun, adaptive.lag = adaptive.lag,
                 fixed.parameters=result$fixed.parameters, 
                 data=result$data)

fin <- try(summary(out), silent=TRUE)

if (inherits(fin, "try-error")) {
  lp <- rep(NA, nrow(out$parametersMCMC$parameters))
  names(lp) <- rownames(out$parametersMCMC$parameters)
  out <- c(out, SD=list(lp))
} else {
  out <- c(out, SD=list(fin$statistics[,"SD"]))
}

out <- addS3Class(out, "mcmcComposite")

return(out)

}
