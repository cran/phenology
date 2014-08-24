#' phenology_MHmcmc runs the Metropolis-Hastings algorithm for data (Bayesian MCMC)
#' @title Run the Metropolis-Hastings algorithm for data
#' @author Marc Girondot
#' @return A list with resultMCMC being mcmc.list object, resultLnL being likelihoods and parametersMCMC being the parameters used
#' @param n.iter Number of iterations for each step
#' @param parametersMCMC A set of parameters used as initial point for searching with information on priors
#' @param result An object obtained after a SearchR fit
#' @param n.chains Number of replicates
#' @param n.adapt Number of iterations before to store outputs
#' @param thin Number of iterations between each stored output
#' @param trace True or False, shows progress
#' @description Run the Metropolis-Hastings algorithm for data.\cr
#' Deeply modified from a MCMC script by Olivier Martin (INRA, Paris-Grignon).\cr
#' The number of iterations is n.iter+n.adapt+1 because the initial likelihood is also displayed.\cr
#' I recommend that thin=1 because the method to estimate SE uses resampling.\cr
#' As initial point is maximum likelihood, n.adapt = 0 seems a good solution.
#' @examples 
#' \dontrun{
#' library(phenology)
#' data(Gratiot)
#' # Generate a formatted list named data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#'     reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' result_Gratiot<-fit_phenology(data=data_Gratiot, 
#' 		parametersfit=parg, parametersfixed=NULL, trace=1)
#' # Generate set of priors for Bayesian analysis
#' pmcmc <- phenology_MHmcmc_p(result_Gratiot, accept = TRUE)
#' result_Gratiot_mcmc <- phenology_MHmcmc(result = result_Gratiot, n.iter = 10000, 
#' parametersMCMC = pmcmc, n.chains = 1, n.adapt = 0, thin = 1, trace = FALSE)
#' # Get standard error of parameters
#' summary(result_Gratiot_mcmc)
#' # Make diagnostics of the mcmc results using coda package
#' mcmc <- as.mcmc(result_Gratiot_mcmc)
#' require(coda)
#' heidel.diag(mcmc)
#' raftery.diag(mcmc)
#' autocorr.diag(mcmc)
#' acf(mcmc[[1]][,"LengthB"], lag.max=200, bty="n", las=1)
#' acf(mcmc[[1]][,"Max_Gratiot"], lag.max=50, bty="n", las=1)
#' batchSE(mcmc, batchSize=100)
#' # The batch standard error procedure is usually thought to 
#' # be not as accurate as the time series methods used in summary
#' summary(mcmc)$statistics[,"Time-series SE"]
#' plot(result_Gratiot_mcmc, parameters=3, las=1, xlim=c(-10, 300))
#' }
#' @export


phenology_MHmcmc<-function(result=stop("An output from fit_phenology() must be provided"), n.iter=10000, 
parametersMCMC=stop("A model generated with phenology_MHmcmc_p() must be provided"), n.chains = 4, 
n.adapt = 0, thin=1, trace=FALSE)
{

  if (class(result)!="phenology") {
    warning("An output from fit_phenology() must be provided")
    return()
  }
  

# likelihood <- .phenology_fonctionMCMC

pt=list(data=result$data, fixed=result$parametersfixed, incertitude=result$method_incertitude, zerocounts=result$zero_counts)


parameters <- parametersMCMC

print(parameters)

out <- .MHalgoGen(n.iter=n.iter, parameters=parametersMCMC, n.chains = n.chains, n.adapt = n.adapt, thin=thin, trace=trace, data=pt, likelihood=.phenology_fonctionMCMC)

return(out)

}
