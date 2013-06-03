\name{likelihood_phenology}
\alias{likelihood_phenology}
\title{Estimate the likelihood of timeseries based on a set of parameters.}
\usage{
  likelihood_phenology(data = NULL, parametersfit = NULL,
    parametersfixed = NULL, zero_counts = NULL,
    method_incertitude = NULL, result = NULL)
}
\arguments{
  \item{data}{Dataset generated with add_format}

  \item{parametersfixed}{Set of fixed parameters}

  \item{parametersfit}{Set of parameters to be fitted}

  \item{method_incertitude}{2 [default] is the correct one
  from a statistical point of view; \cr 0 is an aproximate
  method more rapid; \cr 1 is an alternative more rapid but
  biased.}

  \item{zero_counts}{example c(TRUE, TRUE, FALSE) indicates
  whether the zeros have been recorder for each of these
  timeseries. Defaut is TRUE for all.}

  \item{result}{An object obtained after fit_phenology()}
}
\value{
  The likelihood of the data with the parameters
}
\description{
  This function is used to estimate the likelihood based on
  a set of parameters.
}
\details{
  likelihood_phenology estimate likelihood for a set of
  parameters.
}
\examples{
\dontrun{
# Read a file with data
Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
data(Gratiot)
# Generate a formated list nammed data_Gratiot
data_Gratiot<-add_phenology(Gratiot, name="Complete",
		reference=as.Date("2001-01-01"), format="\%d/\%m/\%Y")
# Generate initial points for the optimisation
parg<-par_init(data_Gratiot, parametersfixed=NULL)
# Estimate likelihood with this initial set of parameters
likelihood_phenology(data=data_Gratiot, parametersfit=parg, parametersfixed=NULL)
}
}
\author{
  Marc Girondot
}

