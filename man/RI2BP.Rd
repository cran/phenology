% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RI2BP.R
\name{RI2BP}
\alias{RI2BP}
\title{Calculate Breeding Proportion from Remigration Interval.}
\usage{
RI2BP(
  proportion = 1,
  RI = c(mean = 2, sd = 0),
  sampling = 10000,
  replicates = 100
)
}
\arguments{
\item{proportion}{A vector of the proportion of different categories}

\item{RI}{The RI of individuals - See description}

\item{sampling}{Number of individuals for sampling}

\item{replicates}{Number of replicates to estimate SE}
}
\value{
Return a vector with mean and se.\cr
}
\description{
This function calculates breeding proportion (BP) from Remigration Interval (RI).\cr
RI can be a vector of RI, one per individual, or an aggregated value with mean and sd and
several categories can exist. See examples.
}
\details{
RI2BP calculate Breeding Proportion from Remigration Interval.
}
\examples{
\dontrun{
library(phenology)
# Example
RI <- c(mean=3, sd=0)
RI2BP(RI=RI)

RI <- c(mean=2, sd=0.3)
RI2BP(RI=RI)

RI <- c(mean=4, sd=0)
RI2BP(RI=RI)

RI <- c(mean=c(2, 10), sd=c(0, 0))
proportion <- c(0.5, 0.5)
RI2BP(proportion=proportion, RI=RI)

c <- c(c1=0.1, c2=0.72, c3=0.126, c4=0.0378, c5=0.0162)
plot(1:5, c, xlab="RI", ylab="Proportion", las=1, bty="n", type="h")
RI2BP(proportion=c, RI=c(mean=1:5, sd=rep(0, 5)))

# To generate random RI with known mean and sd using 
# a truncated lognormal distribution (because 0 does not exist)
mean <- 2.5 - 1
sd <- 1
location <- log(mean^2 / sqrt(sd^2 + mean^2)) 
shape <- sqrt(log(1 + (sd^2 / mean^2)))   
RI <- round(rlnorm(100, meanlog = location, sdlog = shape))+1
mean(RI); sd(RI) # All is ok !
plot(table(RI), xlab="RI", ylab="Proportion", las=1, bty="n", type="h")
RI2BP(proportion=proportion, RI=RI)

}
}
\seealso{
Other Model of Remigration Interval: 
\code{\link{Bayesian.remigration}()},
\code{\link{LnRI_norm}()},
\code{\link{RI}()},
\code{\link{plot.Remigration}()}
}
\author{
Marc Girondot \email{marc.girondot@gmail.com}
}
\concept{Model of Remigration Interval}
