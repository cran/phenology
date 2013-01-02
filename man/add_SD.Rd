\name{add_SD}
\alias{add_SD}
\title{Add SD for a fixed parameter.}
\usage{
  add_SD(parametersfixed = NULL, parameters = NULL,
    SD = NULL, help = FALSE)
}
\arguments{
  \item{parametersfixed}{Set of fixed parameters}

  \item{parameters}{Set of current parameters}

  \item{SD}{Standard deviation value to be added}

  \item{help}{If TRUE, an help is displayed}
}
\value{
  The parameters set with the new SD value
}
\description{
  This function is used to add standard deviation for a
  fixed parameter.
}
\details{
  add_SD adds SD for a fixed parameter.
}
\examples{
library(phenology)
# Generate a set of fixed parameter: Flat and Min
pfixed<-c(Flat=0, Min=0)
	# Add SD for the Flat parameter
pfixed<-add_SD(parametersfixed=pfixed, parameters="Flat", SD=5)
}
\author{
  Marc Girondot
}

