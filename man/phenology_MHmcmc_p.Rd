\name{phenology_MHmcmc_p}
\alias{phenology_MHmcmc_p}
\title{Generates set of parameters to be used with phenology_MHmcmc()}
\usage{
  phenology_MHmcmc_p(result = stop("An output from fit_phenology() must be provided"),
    accept = FALSE)
}
\arguments{
  \item{result}{An object obtained after a SearchR fit}

  \item{accept}{If TRUE, does not wait for}
}
\value{
  A matrix with the parameters
}
\description{
  Interactive script used to generate set of parameters to
  be used with phenology_MHmcmc().\cr
}
\details{
  phenology_MHmcmc_p generates set of parameters to be used
  with MHmcmc()
}
\examples{
library(phenology)
# Read a file with data
\dontrun{
Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
}
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
\dontrun{
pmcmc <- phenology_MHmcmc_p(result_Gratiot)
}
pmcmc <- structure(c("dunif", "dunif", "dunif", "dunif", "dunif", "dunif",
"dunif", "dunif", "0", "0", "0", "0", "0", "0", "0", "0", "200",
"365", "200", "50", "200", "5", "5", "10", "2", "2", "2", "2",
"2", "2", "2", "2", "0", "0", "0", "0", "0", "0", "0", "0", "200",
"365", "200", "50", "200", "5", "5", "10", "95.826796339888",
"175.36499338462", "62.4313052780003", "6.77668901451618e-05",
"33.1138407661406", "0.21779065736816", "0.424368825094697",
"3.58302217559733"), .Dim = c(8L, 7L), .Dimnames = list(c("LengthB",
"Peak", "LengthE", "Flat", "Max_Gratiot", "MinB_Gratiot", "MinE_Gratiot",
"Theta"), c("Density", "Prior1", "Prior2", "SDProp", "Min", "Max",
"Init")))
}
\author{
  Marc Girondot
}

