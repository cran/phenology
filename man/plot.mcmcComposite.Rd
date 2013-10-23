\name{plot.mcmcComposite}
\alias{plot.mcmcComposite}
\title{Plot the result of a MCMC search}
\usage{
  \method{plot}{mcmcComposite} (x, ..., chain = 1,
    parameters = 1, legend = TRUE)
}
\arguments{
  \item{x}{A mcmcComposite object obtained after
  phenology_fonctionMCMC()}

  \item{chain}{The chain to use}

  \item{parameters}{Name of parameters or their number (see
  description)}

  \item{legend}{If FALSE, the legend is not shown}

  \item{...}{Graphical parameters to be send to hist()}
}
\value{
  None
}
\description{
  Plot the result of a MCMC search.\cr The parameters to
  use can be called by:\cr parameters="all"\cr
  parameters=1:4\cr parameters=c("PAR1", "PAR2", "PAR5")\cr
  parameters=c(TRUE, TRUE, FALSE, TRUE)\cr
}
\details{
  plot.mcmcComposite plots the result of a MCMC search
}
\examples{
library(phenology)
# Read a file with data
data(Gratiot)
# Generate a formatted list named data_Gratiot
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
\dontrun{
result_Gratiot_mcmc <- phenology_MHmcmc(result = result_Gratiot, n.iter = 10000,
		parametersMCMC = pmcmc, n.chains = 1, n.adapt = 0, thin = 1, trace = FALSE)
}
data(result_Gratiot_mcmc)
plot(result_Gratiot_mcmc, parameters=3, xlim=c(230, 330))
}
\author{
  Marc Girondot
}

