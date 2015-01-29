#' as.mcmc Extract mcmc object from the result of phenology_MHmcmc to be used with coda package
#' @title Extract mcmc object from the result of phenology_MHmcmc to be used with coda package
#' @author Marc Girondot
#' @return A mcmc.list object
#' @param x A result MHmcmc search
#' @description Take a mcmcComposite object and create a mcmc.list object
#' @examples
#' \dontrun{
#' library(phenology)
#' data(result_Gratiot_mcmc)
#' mcmc <- as.mcmc(result_Gratiot_mcmc)
#' }
#' @import coda
#' @method as.mcmc mcmcComposite
#' @export

as.mcmc.mcmcComposite <-
function(x) {
	if (class(x)!="mcmcComposite") {
		warning("mcmcComposite object must be provided")
	
	} else {

		return(x$resultMCMC)
	}

}
