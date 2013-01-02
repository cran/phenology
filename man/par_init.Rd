\name{par_init}
\alias{par_init}
\title{Calculate initial set of parameters.}
\usage{
  par_init(data = stop("A dataset must be provided"),
    parametersfixed = NA, help = FALSE)
}
\arguments{
  \item{data}{dataset generated with add_format}

  \item{parametersfixed}{Set of fixed parameters}

  \item{help}{If TRUE, an help is displayed}
}
\value{
  The initial set of parameters
}
\description{
  This function is used to generate a first set of
  parameters that is expected to be not to far from the
  final.
}
\details{
  par_init calculates initial set of parameters.
}
\examples{
library(phenology)
# Read a file with data
\dontrun{
Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
}
data(Gratiot)
# Generate a formatted list nammed data_Gratiot
data_Gratiot<-add_phenology(Gratiot, name="Complete", reference=as.Date("2001-01-01"), format="\%d/\%m/\%Y")
# Generate initial points for the optimisation
parg<-par_init(data_Gratiot, parametersfixed=NULL)
# Run the optimisation
\dontrun{
result_Gratiot<-fit_phenology(data=data_Gratiot, parametersfit=parg, parametersfixed=NULL, trace=1)
}
data(result_Gratiot)
# Plot the phenology and get some stats
output<-plot_phenology(result=result_Gratiot, pdf=FALSE)
}
\author{
  Marc Girondot
}

