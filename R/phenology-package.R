#' Fit a parametric function that describes phenology
#'
#' \tabular{ll}{
#'  Package: \tab phenology\cr
#'  Type: \tab Package\cr
#'  Version: \tab 3.54 build 357\cr
#'  Date: \tab 2013-10-23\cr
#'  License: \tab GPL (>= 2)\cr
#'  LazyLoad: \tab yes\cr
#'  }
#' @title The package phenology
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType package
#' @name phenology-package
#' @description A package to fit seasonality counts
#' @references Girondot, M. 2010. Estimating density of animals during 
#'             migratory waves: application to marine turtles at nesting site. 
#'             Endangered Species Research, 12, 85-105.
#' @references Girondot, M. 2010. Editorial: The zero counts. Marine 
#'             Turtle Newsletter, 129, 5-6.
#' @keywords Seasonality Phenology Ecology
#' @seealso Girondot, M., Rivalan, P., Wongsopawiro, R., Briane, J.-P., Hulin, V.,
#'          Caut, S., Guirlet, E. & Godfrey, M. H. 2006. Phenology of marine turtle 
#'          nesting revealed by a statistical model of the nesting season. BMC Ecology, 
#'          6, 11.
#' @examples
#' library(phenology)
#' # Read a file with data
#' \dontrun{
#' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' result_Gratiot<-fit_phenology(data=data_Gratiot, 
#' 		parametersfit=parg, parametersfixed=NULL, trace=1)
#' data(result_Gratiot)
#' # Plot the phenology and get some stats
#' output<-plot(result_Gratiot)
#' }

NULL

