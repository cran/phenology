% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.IP.R
\name{summary.IP}
\alias{summary.IP}
\title{Print the result information from a IP object.}
\usage{
\method{summary}{IP}(object, ..., N = NULL, probs = c(0.025, 0.975))
}
\arguments{
\item{object}{A file of class IP}

\item{...}{Not used}

\item{N}{Number of replicates}

\item{probs}{Probability of confidence interval}
}
\value{
Nothing
}
\description{
The function summary.IP shows result and estimates confidence interval.
}
\details{
summary.IP prints the information from a IP object.
}
\examples{
\dontrun{
library(phenology)
# Read a file with data
}
}
\seealso{
Other Model of Internesting Period: 
\code{\link{IPFit}()},
\code{\link{IPModel}()},
\code{\link{IPPredict}()},
\code{\link{plot.IP}()}
}
\author{
Marc Girondot
}
\concept{Model of Internesting Period}
