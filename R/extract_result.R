#' extract_result get the fitted parameters from a result object.
#' @title Extract the set of parameters from a result object.
#' @author Marc Girondot
#' @return Return the set of fitted parameters
#' @param result A result file
#' @description The function "extract_result" permits to extract the set of parameters from a result object obtained after fit_phenology.
#' @family Phenology model
#' @examples
#' library(phenology)
#' \dontrun{
#' # Read a file with data
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, fixed.parameters=NULL)
#' # Run the optimisation
#' # result_Gratiot<-fit_phenology(data=data_Gratiot, fitted.parameters=parg, 
#' 		fixed.parameters=NULL)
#' data(result_Gratiot)
#' # Extract the fitted parameters
#' parg1<-extract_result(result_Gratiot)
#' }
#' @export


extract_result <-
function(result=NULL) {
	return(result$par)
}
