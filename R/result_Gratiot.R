#' Result of the fit of Leatherback nest counts from Gratiot et al. (2006) Figure 1
#' @title Result of the fit of Leatherback nest counts
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name result_Gratiot
#' @description Result of the fit of Leatherback nest counts from Gratiot et al. (2006) Figure 1
#' The phenology has been fitted with MinE, MinB, Max, Flat, LengthB, LengthE, Peak, Theta.
#' @references Gratiot, N., Gratiot, J., de Thoisy, B. & Kelle, L. 2006. 
#'             Estimation of marine turtles nesting season from incomplete 
#'             data ; statistical adjustment of a sinusoidal function. Animal 
#'             Conservation, 9, 95-102.
#' @keywords datasets
#' @usage result_Gratiot
#' @examples
#' \dontrun{
#' library(phenology)
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot <- add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg <- par_init(data_Gratiot, fixed.parameters=NULL)
#' # Run the optimisation
#' result_Gratiot <- fit_phenology(data=data_Gratiot, 
#' 		fitted.parameters=parg, fixed.parameters=NULL)
#' # Read a file with result
#' data(result_Gratiot)
#' }
#' @format A list with Gratiot data and the result of the fit.
NULL
