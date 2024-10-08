% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Tagloss_LengthObs.R
\name{Tagloss_LengthObs}
\alias{Tagloss_LengthObs}
\title{Return a list with the number of days for different kinds of individuals are seen.}
\usage{
Tagloss_LengthObs(data, progressbar = TRUE)
}
\arguments{
\item{data}{Set of indivuals}

\item{progressbar}{Is shown a progressbar?}
}
\value{
Return a list with the number of days for different kinds of individuals are seen.
}
\description{
Usefull to summarize data
}
\details{
Tagloss_LengthObs returns a list with the number of days for different kinds of individuals are seen.
}
\examples{
\dontrun{
library(phenology)
# Example
data_f_21 <- Tagloss_format(outLR, model="21")
Tagloss_LengthObs(data_f_21)
}
}
\seealso{
Other Model of Tag-loss: 
\code{\link{Tagloss_L}()},
\code{\link{Tagloss_cumul}()},
\code{\link{Tagloss_daymax}()},
\code{\link{Tagloss_fit}()},
\code{\link{Tagloss_format}()},
\code{\link{Tagloss_mcmc}()},
\code{\link{Tagloss_mcmc_p}()},
\code{\link{Tagloss_model}()},
\code{\link{Tagloss_simulate}()},
\code{\link{logLik.Tagloss}()},
\code{\link{o_4p_p1p2}},
\code{\link{plot.Tagloss}()},
\code{\link{plot.TaglossData}()}
}
\author{
Marc Girondot \email{marc.girondot@gmail.com}
}
\concept{Model of Tag-loss}
