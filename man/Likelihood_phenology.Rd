\name{likelihood_phenology}
\alias{likelihood_phenology}
\title{Estimate the likelihood of timeseries based on a set of parameters.}
\usage{
  likelihood_phenology(data = dta, parametersfit = x,
    parametersfixed = NULL, zero_counts = TRUE,
    method_incertitude = 0, help = FALSE)
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

  \item{help}{If TRUE, an help is displayed}
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
# Read a file with data
# Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", , header=FALSE)
data(Gratiot)
# Generate a formatted list nammed data_Gratiot
data_Gratiot<-add_format(origin=NULL, add=Gratiot, name="Complete", reference=as.Date("2001-01-01"), format="\%d/\%m/\%Y")
# Generate initial points for the optimisation
parg<-par_init(data_Gratiot, parametersfixed=NULL)
# Estimate likelihood with this initial set of parameters
likelihood_phenology(data=data_Gratiot, parametersfit=parg, parametersfixed=NULL)
}
\author{
  Marc Girondot
}

