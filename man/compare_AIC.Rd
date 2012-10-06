\name{compare_AIC}
\alias{compare_AIC}
\title{Compares the AIC of several outputs}
\usage{
  compare_AIC(result = NULL, help = FALSE)
}
\arguments{
  \item{result}{A list with the results to be compared}

  \item{help}{If TRUE, an help is displayed}
}
\value{
  A list with DeltaAIC and Akaike weight for the models.
}
\description{
  This function is used to compares the AIC of several
  outputs obtained with the same data but with different
  set of parameters.
}
\details{
  compare_AIC compares the AIC of several outputs obtained
  with the same data.
}
\examples{
# Read a file with data
library("phenology")
#' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", , header=FALSE)
data(Gratiot)
# Generate a formated list nammed data_Gratiot
data_Gratiot<-add_format(origin=NULL, add=Gratiot, name="Complete", reference=as.Date("2001-01-01"), format="\%d/\%m/\%Y")
# Fix parameter FLat to 0
pfixed=c(Flat=0)
# Generate initial points for the optimisation
parg<-par_init(data_Gratiot, parametersfixed=pfixed)
# Fit is done
# result_Gratiot_Flat<-fit_phenology(data=data_Gratiot, parametersfit=parg2, parametersfixed=pfixed, trace=1)
data(result_Gratiot_Flat)
# Generate initial points for the optimisation
parg<-par_init(data_Gratiot, parametersfixed=NULL)
# Run the optimisation
# result_Gratiot<-fit_phenology(data=data_Gratiot, parametersfit=parg, parametersfixed=NULL, trace=1)
data(result_Gratiot)
# Compare both models
outputAIC<-compare_AIC(list(full=result_Gratiot, Flat=result_Gratiot_Flat))
}
\author{
  Marc Girondot
}

