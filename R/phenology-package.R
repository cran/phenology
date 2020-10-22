#' Fit a parametric function that describes phenology
#'
#' \tabular{ll}{
#'  Package: \tab phenology\cr
#'  Type: \tab Package\cr
#'  Version: \tab 7.6 build 1168\cr
#'  Date: \tab 2020-10-22\cr
#'  License: \tab GPL (>= 2)\cr
#'  LazyLoad: \tab yes\cr
#'  }
#' @title Tools to Manage a Parametric Function that Describes Phenology and More
#' @author Marc Girondot \email{marc.girondot@@universite-paris-saclay.fr}
#' @docType package
#' @name phenology-package
#' @description Functions used to fit and test the phenology of species based on counts.\cr
#' Note that only the most significant changes are reported in the NEWS.\cr
#' The lastest version of this package can always been installed using:\cr
#' install.packages("http://www.ese.u-psud.fr/epc/conservation/CRAN/HelpersMG.tar.gz", repos=NULL, type="source")\cr
#' install.packages("http://www.ese.u-psud.fr/epc/conservation/CRAN/phenology.tar.gz", repos=NULL, type="source")
#' @references Girondot, M. 2010. Estimating density of animals during 
#'             migratory waves: application to marine turtles at nesting site. 
#'             Endangered Species Research, 12, 85-105.
#' @references Girondot M. and Rizzo A. 2015. Bayesian framework to integrate 
#'             traditional ecological knowledge into ecological modeling: A case 
#'             study. Journal of Ethnobiology, 35, 339-355.
#' @references Girondot, M. 2010. Editorial: The zero counts. Marine 
#'             Turtle Newsletter, 129, 5-6.
#' @references Girondot, M., 2017. Optimizing sampling design to infer marine turtles 
#'             seasonal nest number for low-and high-density nesting beach using convolution 
#'             of negative binomial distribution. Ecological Indicators 81, 83–89.
#' @references Rivalan, P., Godfrey, M.H., Prévot-Julliard, A.-C., Girondot, M., 2005. 
#'             Maximum likelihood estimates of tag loss in leatherback sea turtles. Journal 
#'             of Wildlife Management 69, 540-548.
#' @keywords Seasonality Phenology Ecology tagloss OCF ECF Clutch
#' @seealso Girondot, M., Rivalan, P., Wongsopawiro, R., Briane, J.-P., Hulin, V.,
#'          Caut, S., Guirlet, E. & Godfrey, M. H. 2006. Phenology of marine turtle 
#'          nesting revealed by a statistical model of the nesting season. BMC Ecology, 
#'          6, 11.
#' @seealso Delcroix, E., Bédel, S., Santelli, G., Girondot, M., 2013. Monitoring 
#'          design for quantification of marine turtle nesting with limited human 
#'          effort: a test case in the Guadeloupe Archipelago. Oryx 48, 95-105.
#' @seealso Briane J-P, Rivalan P, Girondot M (2007) The inverse problem applied 
#'             to the Observed Clutch Frequency of Leatherbacks from Yalimapo beach, 
#'             French Guiana. Chelonian Conservation and Biology 6:63-69
#' @seealso Fossette S, Kelle L, Girondot M, Goverse E, Hilterman ML, Verhage B, 
#'          Thoisy B, de, Georges J-Y (2008) The world's largest leatherback 
#'          rookeries: A review of conservation-oriented research in French 
#'          Guiana/Suriname and Gabon. Journal of Experimental Marine Biology 
#'          and Ecology 356:69-82
#' @examples
#' \dontrun{
#' library(phenology)
#' # Get the lastest version at:
#' # install.packages("http://www.ese.u-psud.fr/epc/conservation/CRAN/phenology.tar.gz", 
#' #     repos=NULL, type="source")
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
#' data(result_Gratiot)
#' # Plot the phenology and get some stats
#' output <- plot(result_Gratiot)
#' }

NULL

