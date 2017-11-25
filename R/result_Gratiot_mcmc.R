#' Result of the mcmc for Leatherback nest counts from Gratiot et al. (2006) Figure 1
#' @title Result of the mcmc for Leatherback nest counts
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name result_Gratiot_mcmc
#' @description Result of the mcmc for Leatherback nest counts from Gratiot et al. (2006) Figure 1
#' The phenology has been fitted with MinE, MinB, Max, Flat, LengthB, LengthE, Peak, Theta.
#' @references Gratiot, N., Gratiot, J., de Thoisy, B. & Kelle, L. 2006. 
#'             Estimation of marine turtles nesting season from incomplete 
#'             data ; statistical adjustment of a sinusoidal function. Animal 
#'             Conservation, 9, 95-102.
#' @keywords datasets
#' @usage result_Gratiot_mcmc
#' @examples
#' \dontrun{
#' library(phenology)
#' data(result_Gratiot)
#' # Read a file with data
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot <- add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg <- par_init(data_Gratiot, fixed.parameters=NULL)
#' # Run the optimisation
#' result_Gratiot <- fit_phenology(data=data_Gratiot, 
#' 		fitted.parameters=parg, fixed.parameters=NULL)
#' # generate data for mcmc run
#' pmcmc <- phenology_MHmcmc_p(result_Gratiot, accept = TRUE)
#' result_Gratiot_mcmc <- phenology_MHmcmc(result = result_Gratiot, 
#'      n.iter = 10000, 
#'      adaptive=TRUE,
#'      parametersMCMC = pmcmc, 
#'      n.chains = 1, n.adapt = 0, thin = 1, trace = FALSE)
#' # Read a file with result
#' data(result_Gratiot_mcmc)
#' 1-rejectionRate(as.mcmc(result_Gratiot_mcmc))
#' 
#' summary(result_Gratiot, resultmcmc=result_Gratiot_mcmc)
#' }
#' @format A mcmcComposite object with mcmc result.
NULL
