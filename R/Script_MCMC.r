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
#' library(phenology)
#' # Read a file with data
#' # Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", , header=FALSE)
#' data(Gratiot)
#' # Generate a formatted list named data_Gratiot 
#' data_Gratiot<-add_format(origin=NULL, add=Gratiot, name="Complete", reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' ## not run
#' # result_Gratiot<-fit_phenology(data=data_Gratiot, parametersfit=parg, parametersfixed=NULL, trace=1)
#' ## end not run
#' data(result_Gratiot)
#' # Generate set of priors for Bayesian analysis
#' ## not run
#' ## pmcmc <- phenology_MHmcmc_p(result_Gratiot)
#' ## end not run
#' pmcmc <- structure(c("dunif", "dunif", "dunif", "dunif", "dunif", "dunif", 
#' "dunif", "dunif", "0", "0", "0", "0", "0", "0", "0", "0", "200", 
#' "365", "200", "50", "200", "5", "5", "10", "2", "2", "2", "2", 
#' "2", "2", "2", "2", "0", "0", "0", "0", "0", "0", "0", "0", "200", 
#' "365", "200", "50", "200", "5", "5", "10", "95.826796339888", 
#' "175.36499338462", "62.4313052780003", "6.77668901451618e-05", 
#' "33.1138407661406", "0.21779065736816", "0.424368825094697", 
#' "3.58302217559733"), .Dim = c(8L, 7L), .Dimnames = list(c("LengthB", 
#' "Peak", "LengthE", "Flat", "Max_Gratiot", "MinB_Gratiot", "MinE_Gratiot", 
#' "Theta"), c("Density", "Prior1", "Prior2", "SDProp", "Min", "Max", 
#' "Init")))
#' ## not run
#' # res_mcmc <- phenology_MHmcmc(result = result_Gratiot, n.iter = 10000, 
#' # parametersMCMC = pmcmc, n.chains = 1, n.adapt = 0, thin = 1, trace = FALSE)
#' ## end not run
#' @export


phenology_MHmcmc<-function(result=stop("An output from searchR must be provided"), n.iter=10000, 
parametersMCMC=stop("A model generated with phenology_MHmcmc_p() must be provided"), n.chains = 4, n.adapt = 0, thin=1, trace=FALSE)
{


# likelihood <- .phenology_fonctionMCMC

pt=list(data=result$data, fixed=result$parametersfixed, incertitude=result$method_incertitude, zerocounts=result$zero_counts)


parameters <- parametersMCMC

print(parameters)

out <- .MHalgoGen(n.iter=n.iter, parameters=parametersMCMC, n.chains = n.chains, n.adapt = n.adapt, thin=thin, trace=trace, data=pt, likelihood=.phenology_fonctionMCMC)

return(out)

}
