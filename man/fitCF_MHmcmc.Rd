% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitCF_MHmcmc.R
\name{fitCF_MHmcmc}
\alias{fitCF_MHmcmc}
\title{Run the Metropolis-Hastings algorithm for ECFOCF data}
\usage{
fitCF_MHmcmc(
  result = stop("An output from fitCF() must be provided"),
  n.iter = 10000,
  parametersMCMC = stop("A parameter set from fitCF_MHmcmc_p() must be provided"),
  n.chains = 1,
  n.adapt = 0,
  thin = 1,
  adaptive = FALSE,
  adaptive.lag = 500,
  adaptive.fun = function(x) {
     ifelse(x > 0.234, 1.3, 0.7)
 },
  trace = FALSE,
  traceML = FALSE,
  intermediate = NULL,
  filename = "intermediate.Rdata",
  previous = NULL
)
}
\arguments{
\item{result}{An object obtained after a SearchR fit}

\item{n.iter}{Number of iterations for each step}

\item{parametersMCMC}{A set of parameters used as initial point for searching with information on priors}

\item{n.chains}{Number of replicates}

\item{n.adapt}{Number of iterations before to store outputs}

\item{thin}{Number of iterations between each stored output}

\item{adaptive}{Should an adaptive process for SDProp be used}

\item{adaptive.lag}{Lag to analyze the SDProp value in an adaptive content}

\item{adaptive.fun}{Function used to change the SDProp}

\item{trace}{TRUE or FALSE or period, shows progress}

\item{traceML}{TRUE or FALSE to show ML}

\item{intermediate}{Period for saving intermediate result, NULL for no save}

\item{filename}{If intermediate is not NULL, save intermediate result in this file}

\item{previous}{Previous result to be continued. Can be the filename in which intermediate results are saved.}
}
\value{
A list with resultMCMC being mcmc.list object, resultLnL being likelihoods and parametersMCMC being the parameters used
}
\description{
Run the Metropolis-Hastings algorithm for RMU.data.\cr
The number of iterations is n.iter+n.adapt+1 because the initial likelihood is also displayed.\cr
I recommend thin=1 because the method to estimate SE uses resampling.\cr
As initial point is maximum likelihood, n.adapt = 0 is a good solution.\cr
The parameters intermediate and filename are used to save intermediate results every 'intermediate' iterations (for example 1000). Results are saved in a file of name filename.\cr
The parameter previous is used to indicate the list that has been save using the parameters intermediate and filename. It permits to continue a mcmc search.\cr
These options are used to prevent the consequences of computer crash or if the run is very very long and computer processes at time limited.\cr
}
\details{
fitCF_MHmcmc runs the Metropolis-Hastings algorithm for ECFOCF (Bayesian MCMC)
}
\examples{
\dontrun{
library("phenology")
data(MarineTurtles_2002)
ECFOCF_2002 <- TableECFOCF(MarineTurtles_2002)

# Paraetric model for clutch frequency
o_mu1p1_CFp <- fitCF(x = c(mu = 2.1653229641404539, 
                 sd = 1.1465246643327098, 
                 p = 0.25785366120357966), 
                 fixed.parameters=NULL, 
                 data=ECFOCF_2002, hessian = TRUE)
                           
pMCMC <- fitCF_MHmcmc_p(result=o_mu1p1_CFp, accept=TRUE)
fitCF_MCMC <- fitCF_MHmcmc(result = o_mu1p1_CFp, n.iter = 1000, 
                           parametersMCMC = pMCMC, n.chains = 1, n.adapt = 0, 
                           adaptive=TRUE, 
                           thin = 1, trace = TRUE)
                           
plot(fitCF_MCMC, parameters="mu")
plot(fitCF_MCMC, parameters="sd")
plot(fitCF_MCMC, parameters="p", xlim=c(0, 0.5), breaks=seq(from=0, to=0.5, by=0.05))
plot(fitCF_MCMC, parameters="p", transform = invlogit, xlim=c(0, 1), 
     breaks=c(seq(from=0, to=1, by=0.05)))

}
}
\seealso{
Other Model of Clutch Frequency: 
\code{\link{ECFOCF_f}()},
\code{\link{ECFOCF_full}()},
\code{\link{TableECFOCF}()},
\code{\link{fitCF}()},
\code{\link{fitCF_MHmcmc_p}()},
\code{\link{generateCF}()},
\code{\link{lnLCF}()},
\code{\link{logLik.ECFOCF}()},
\code{\link{plot.ECFOCF}()},
\code{\link{plot.TableECFOCF}()}
}
\author{
Marc Girondot \email{marc.girondot@gmail.com}
}
\concept{Model of Clutch Frequency}
