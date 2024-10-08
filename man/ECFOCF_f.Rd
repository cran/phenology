% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ECFOCF_f.R
\name{ECFOCF_f}
\alias{ECFOCF_f}
\title{Calculate a table of probabilities of ECF and OCF.}
\usage{
ECFOCF_f(
  mu,
  sd = NA,
  p,
  MaxNests = 15,
  mu_season = NA,
  sd_season = NA,
  MeanDaysBetween2Nests = 9.8,
  length_season = floor(365/MeanDaysBetween2Nests) + 1,
  parallel = TRUE
)
}
\arguments{
\item{mu}{The average of lognormal for clutch frequency.}

\item{sd}{The sd parameter of lognormal for clutch frequency.}

\item{p}{The capture probability for an individual nesting event. As a probability.}

\item{MaxNests}{Maximum number of nests by a female.}

\item{mu_season}{The average of ordinal day for beginning of nesting season.}

\item{sd_season}{The sd parameter of lognormal for ordinal day for beginning of nesting season.}

\item{MeanDaysBetween2Nests}{Average number of days between two nests.}

\item{length_season}{The total length of season based on groups of interclutch intervals.}

\item{parallel}{If TRUE parallel computing is used.}
}
\value{
Return a matrix of class TableECFOCF.\cr
}
\description{
This function calculates a table of probabilities of ECF and OCF.\cr
If p is lower or higher than 1E-100 or 1-1E-100, it is changed to 1E-100 and 1-(1E-100) respectively.\cr
Names for p vector elements should be p, or px (with x=1:categories), or px.period.\cr
If mu_season and sd_season are equal to NA, the model is not temporalized.\cr
If mu_season and sd_season are not NA, the model returns a 3D-table OCFECF.\cr
}
\details{
ECFOCF_f calculate a table of probabilities of ECF and OCF.
}
\examples{
\dontrun{
library(phenology)
# Example
modelECFOCF <- ECFOCF_f(mu=5.58013243236187, 
                    sd=1.225581130238, 
                    p=invlogit(1.3578137414575), 
                    MaxNests=15)
plot(modelECFOCF)
modelECFOCF <- ECFOCF_f(mu=5.58013243236187, 
                    sd=1.225581130238, 
                    mu_season=12, 
                    sd_season=2, 
                    p=c(p1=invlogit(1.3578137414575)), 
                    MaxNests=15, 
                    MeanDaysBetween2Nests=9.8, 
                    length_season=floor(365/9.8)+1)
plot(modelECFOCF, period=2)
}
}
\seealso{
Other Model of Clutch Frequency: 
\code{\link{ECFOCF_full}()},
\code{\link{TableECFOCF}()},
\code{\link{fitCF}()},
\code{\link{fitCF_MHmcmc}()},
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
