\name{as.mcmc.mcmcComposite}
\alias{as.mcmc.mcmcComposite}
\title{Extract mcmc object from the result of MHmcmc to be used with coda package}
\usage{
\method{as.mcmc}{mcmcComposite}(x)
}
\arguments{
  \item{x}{A result MHmcmc search}
}
\value{
A mcmc.list object
}
\description{
Take a mcmcComposite object and create a mcmc.list object
}
\details{
as.mcmc Extract mcmc object from the result of MHmcmc to be
used with coda package
}
\examples{
\dontrun{
library(phenology)
data(result_Gratiot_mcmc)
mcmc <- as.mcmc(result_Gratiot_mcmc)
}
}
\author{
Marc Girondot
}

