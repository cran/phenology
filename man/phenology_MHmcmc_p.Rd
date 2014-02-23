\name{phenology_MHmcmc_p}
\alias{phenology_MHmcmc_p}
\title{Generates set of parameters to be used with phenology_MHmcmc()}
\usage{

  phenology_MHmcmc_p(result = stop("An output from fit_phenology() must be provided"),
  accept = FALSE)
}
\arguments{
  \item{result}{An object obtained after a fit_phenology()
  fit}

  \item{accept}{If TRUE, does not wait for use interaction}
}
\value{
A matrix with the parameters
}
\description{
Interactive script used to generate set of parameters to be
used with phenology_MHmcmc().\cr
}
\details{
phenology_MHmcmc_p generates set of parameters to be used
with MHmcmc()
}
\examples{
library(phenology)
# Read a file with data
data(Gratiot)
# Generate a formatted list nammed data_Gratiot
data_Gratiot<-add_phenology(Gratiot, name="Complete",
		reference=as.Date("2001-01-01"), format="\%d/\%m/\%Y")
# Generate initial points for the optimisation
parg<-par_init(data_Gratiot, parametersfixed=NULL)
# Run the optimisation
\dontrun{
result_Gratiot<-fit_phenology(data=data_Gratiot,
		parametersfit=parg, parametersfixed=NULL, trace=1)
}
data(result_Gratiot)
# Generate set of priors for Bayesian analysis
pmcmc <- phenology_MHmcmc_p(result_Gratiot, accept=TRUE)
}
\author{
Marc Girondot
}

