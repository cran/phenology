% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitRMU_MHmcmc_p.R
\name{fitRMU_MHmcmc_p}
\alias{fitRMU_MHmcmc_p}
\title{Generates set of parameters to be used with fitRMU_MHmcmc()}
\usage{
fitRMU_MHmcmc_p(
  result = stop("An output from fitRMU() must be provided"),
  density = "dunif",
  accept = FALSE
)
}
\arguments{
\item{result}{An object obtained after a fitRMU() fit}

\item{density}{Preset of density; can be dnorm, dunif, or dgamma}

\item{accept}{If TRUE, does not wait for user interaction}
}
\value{
A matrix with the parameters
}
\description{
Interactive or automatic script used to generate set of parameters to be
used with fitRMU_MHmcmc().\cr
If density="dgamma" is used, a uniform distribution is used for r, as
r can be negative.
}
\details{
fitRMU_MHmcmc_p generates set of parameters to be used with fitRMU_MHmcmc()
}
\examples{
\dontrun{
library("phenology")
RMU.name.AtlanticW <- data.frame(mean=c("Yalimapo.French.Guiana", 
                                         "Galibi.Suriname", 
                                         "Irakumpapy.French.Guiana"), 
                                 se=c("se_Yalimapo.French.Guiana", 
                                      "se_Galibi.Suriname", 
                                      "se_Irakumpapy.French.Guiana"))
data.AtlanticW <- data.frame(Year=c(1990:2000), 
      Yalimapo.French.Guiana=c(2076, 2765, 2890, 2678, NA, 
                               6542, 5678, 1243, NA, 1566, 1566),
      se_Yalimapo.French.Guiana=c(123.2, 27.7, 62.5, 126, NA, 
                                 230, 129, 167, NA, 145, 20),
      Galibi.Suriname=c(276, 275, 290, NA, 267, 
                       542, 678, NA, 243, 156, 123),
      se_Galibi.Suriname=c(22.3, 34.2, 23.2, NA, 23.2, 
                           4.3, 2.3, NA, 10.3, 10.1, 8.9),
      Irakumpapy.French.Guiana=c(1076, 1765, 1390, 1678, NA, 
                               3542, 2678, 243, NA, 566, 566),
      se_Irakumpapy.French.Guiana=c(23.2, 29.7, 22.5, 226, NA, 
                                 130, 29, 67, NA, 15, 20))
                           
cst <- fitRMU(data=data.AtlanticW, RMU.name=RMU.name.AtlanticW, 
               colname.year="Year", model.trend="Constant", 
               model.SD="Zero")
pMCMC <- fitRMU_MHmcmc_p(result=cst, accept=TRUE)
}
}
\seealso{
Other Fill gaps in RMU: 
\code{\link{CI.RMU}()},
\code{\link{fitRMU}()},
\code{\link{fitRMU_MHmcmc}()},
\code{\link{logLik.fitRMU}()},
\code{\link{plot.fitRMU}()}
}
\author{
Marc Girondot
}
\concept{Fill gaps in RMU}
