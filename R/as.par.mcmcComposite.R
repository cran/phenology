#' as.par.mcmcComposite Extract parameters at maximum likelihood from the result of MHmcmc
#' @title Extract parameters at maximum likelihood from the result of MHmcmc
#' @author Marc Girondot
#' @return A vector with parameters at maximum likelihood
#' @param x A result from MHmcmc search
#' @description Take a mcmcComposite object and create a vector object with parameter value at maximum likelihood.
#' It also says at which iteration the maximum lihelihood has been observed.
#' @examples
#' \dontrun{
#' library(phenology)
#' data(result_Gratiot_mcmc)
#' x <- as.par.mcmcComposite(result_Gratiot_mcmc)
#' }
#' @export


as.par.mcmcComposite <-
function(x=stop("A result obtained from MHmcmc must be provided")) {

	L <- x$resultLnL[[1]]
	p <- x$resultMCMC[,][[1]]
	
	tab <- 0

	pos <- which.max(L)
	print(paste("The best likelihood has been observed at iteration", pos-tab))
	pml <- as.numeric(p[pos,])
	names(pml) <- names(p[pos,])
	return(pml)

}
