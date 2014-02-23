\name{adapt_parameters}
\alias{adapt_parameters}
\title{Extract the parameters from a set of parameters to be used with another dataset.}
\usage{
adapt_parameters(data = stop("Datasets is mandatory for this function"),
  parameters = stop("Set of parameters is mandatory for this function"))
}
\arguments{
  \item{parameters}{A set of parameters}

  \item{data}{A dataset of counts}
}
\value{
Return the set of parameters
}
\description{
The function "adapt_parameters" extracts the set of
parameters to be used with a subset of data. All the
uncessary parameters are removed. It can be used when a set
of beaches are fitted first and after only one of these
beaches is fitted again.
}
\details{
adapt_parameters get the fitted parameters from a result
object.
}
\examples{
library(phenology)
# Read a file with data
\dontrun{
Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
}
data(Gratiot)
# Generate a formatted list nammed data_Gratiot
refdate <- as.Date("2001-01-01")
data_Gratiot<-add_phenology(Gratiot, name="Complete",
		reference=refdate, format="\%d/\%m/\%Y")
# Generate initial points for the optimisation
parg<-par_init(data_Gratiot, parametersfixed=NULL)
# Add unnecessary parameters to parg
parg <- c(parg, Max_dummybeach=2, Peak_dummybeach=123)
# Extract the fitted parameters
parg1<-adapt_parameters(data=data_Gratiot, parameters=parg)
}
\author{
Marc Girondot
}

