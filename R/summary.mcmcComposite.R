#' summary.mcmcComposite get info on the result of a MCMC search
#' @title Summarize the result of a MCMC search
#' @author Marc Girondot
#' @return A summary of the result
#' @param object A mcmcComposite object obtained after MHmcmc()
#' @param ... Internal use
#' @param chain The chain to use
#' @description Summary for the result of a MCMC search
#' @examples 
#' library(phenology)
#' # Read a file with data
#' \dontrun{
#' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
#' }
#' data(Gratiot)
#' # Generate a formatted list named data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' \dontrun{
#' result_Gratiot<-fit_phenology(data=data_Gratiot, 
#' 		parametersfit=parg, parametersfixed=NULL, trace=1)
#' }
#' data(result_Gratiot)
#' # Generate set of priors for Bayesian analysis
#' \dontrun{
#' pmcmc <- phenology_MHmcmc_p(result_Gratiot)
#' }
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
#' \dontrun{
#' res_mcmc <- phenology_MHmcmc(result = result_Gratiot, n.iter = 10000, 
#' parametersMCMC = pmcmc, n.chains = 1, n.adapt = 0, thin = 1, trace = FALSE)
#' data(res_mcmc)
#' summary(res_mcmc)
#' }
#' @method summary mcmcComposite
#' @export

summary.mcmcComposite <- function(object, ... , chain=NULL) {

resultMCMC <- object

mcmc <- resultMCMC[["resultMCMC"]]

if (!is.null(chain)) {mcmc <- mcmc[[chain]]}

return(summary(mcmc))


}
