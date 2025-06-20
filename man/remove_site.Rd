% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_site.R
\name{remove_site}
\alias{remove_site}
\title{Removes site information from a set of parameters.}
\usage{
remove_site(parameters = NULL, help = FALSE)
}
\arguments{
\item{parameters}{Set of parameters}

\item{help}{If TRUE, an help is displayed}
}
\value{
Return a set of modified parameters
}
\description{
This function is used to remove the information of the site
from a set of parameters. It can be used to other timeseries after.
}
\details{
remove_site removes beach information from a set of parameters.
}
\examples{
library(phenology)
# Read a file with data
data(Gratiot)
# Generate a formatted list nammed data_Gratiot 
data_Gratiot<-add_phenology(Gratiot, name="Complete", 
		reference=as.Date("2001-01-01"), format="\%d/\%m/\%Y")
# Generate initial points for the optimisation
parg<-par_init(data_Gratiot, fixed.parameters=NULL)
# Run the optimisation
\dontrun{
result_Gratiot<-fit_phenology(data=data_Gratiot, 
		fitted.parameters=parg, fixed.parameters=NULL)
}
data(result_Gratiot)
# Extract parameters form result
parg<-extract_result(result_Gratiot)
# Remove site information
parg1<-remove_site(parg)
}
\seealso{
Other Phenology model: 
\code{\link{AutoFitPhenology}()},
\code{\link{BE_to_LBLE}()},
\code{\link{Gratiot}},
\code{\link{LBLE_to_BE}()},
\code{\link{LBLE_to_L}()},
\code{\link{L_to_LBLE}()},
\code{\link{MarineTurtles_2002}},
\code{\link{MinBMinE_to_Min}()},
\code{\link{adapt_parameters}()},
\code{\link{add_SE}()},
\code{\link{add_phenology}()},
\code{\link{extract_result}()},
\code{\link{fit_phenology}()},
\code{\link{likelihood_phenology}()},
\code{\link{logLik.phenology}()},
\code{\link{map_Gratiot}},
\code{\link{map_phenology}()},
\code{\link{par_init}()},
\code{\link{phenology}()},
\code{\link{phenology2fitRMU}()},
\code{\link{phenology_MHmcmc}()},
\code{\link{phenology_MHmcmc_p}()},
\code{\link{plot.phenology}()},
\code{\link{plot.phenologymap}()},
\code{\link{plot_delta}()},
\code{\link{plot_phi}()},
\code{\link{print.phenology}()},
\code{\link{print.phenologymap}()},
\code{\link{print.phenologyout}()},
\code{\link{result_Gratiot}},
\code{\link{result_Gratiot1}},
\code{\link{result_Gratiot2}},
\code{\link{result_Gratiot_Flat}},
\code{\link{summary.phenology}()},
\code{\link{summary.phenologymap}()},
\code{\link{summary.phenologyout}()}
}
\author{
Marc Girondot
}
\concept{Phenology model}
