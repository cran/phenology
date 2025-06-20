% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AutoFitPhenology.R
\name{AutoFitPhenology}
\alias{AutoFitPhenology}
\title{Automatic fit for phenology and tests}
\usage{
AutoFitPhenology(
  data = stop("A dataset must be provided"),
  progressbar = TRUE,
  ...
)
}
\arguments{
\item{data}{Dataset generated with add_phenology()}

\item{progressbar}{If FALSE, do not show the progress bar}

\item{...}{Parameters for fit_phenology()}
}
\value{
A list with 12 elements corresponding to the 12 tested models
}
\description{
This function is used to test several combinations of fit at a time.
}
\details{
AutoFitPhenology runs fit for phenology and tests several combinations
}
\examples{
\dontrun{
library(phenology)
# Read a file with data
data(Gratiot)
# Generate a formatted list nammed data_Gratiot 
data_Gratiot <- add_phenology(Gratiot, name="Complete", 
		reference=as.Date("2001-01-01"), format="\%d/\%m/\%Y")
# Run the optimisation
result_Gratiot_Auto <- AutoFitPhenology(data=data_Gratiot)
result_Gratiot_Auto <- AutoFitPhenology(data=data_Gratiot, 
          control=list(trace=0, REPORT=100, maxit=500))
}
}
\seealso{
Other Phenology model: 
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
\code{\link{remove_site}()},
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
