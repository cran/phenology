\name{map_phenology}
\alias{map_phenology}
\title{Generate a likelihood map varying Phi and Delta.}
\usage{
  map_phenology(data = NULL, parametersfit = NULL,
    parametersfixed = NA,
    Phi = seq(from = 0.2, to = 20, length.out = 100),
    Delta = NULL, method_incertitude = 2,
    zero_counts = TRUE, progressbar = TRUE, help = FALSE)
}
\arguments{
  \item{data}{dataset generated with add_format}

  \item{parametersfixed}{Set of fixed parameters}

  \item{parametersfit}{Set of parameters to be fitted}

  \item{Phi}{Phi values to be analyzed}

  \item{Delta}{Delta value to be analyzed}

  \item{method_incertitude}{2 [default] is the correct one
  from a statistical point of view; \cr 0 is an aproximate
  method more rapid; \cr 1 is an alternative more rapid but
  biased.\cr}

  \item{zero_counts}{Example c(TRUE, TRUE, FALSE) indicates
  whether the zeros have been recorder for each of these
  timeseries. Defaut is TRUE for all.}

  \item{progressbar}{If FALSE, do not show the progress
  bar}

  \item{help}{If TRUE, an help is displayed}
}
\value{
  Display a likelihood map
}
\description{
  This function generates a map of likelihood varying Phi
  and Delta.\cr When Delta is not given, the same precision
  as Phi is used.
}
\details{
  map_phenology generates a likelihood map.
}
\examples{
library("phenology")
# Read a file with data
\dontrun{
Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
}
data(Gratiot)
# Generate a formatted list nammed data_Gratiot
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
# Extract the fitted parameters
parg1<-extract_result(result_Gratiot)
# Add constant Alpha and Tau values
# [day d amplitude=(Alpha+Nd*Beta)^Tau with Nd being the number of counts for day d]
pfixed<-c(parg1, Alpha=0, Tau=1)
pfixed<-pfixed[-which(names(pfixed)=="Theta")]
# The only fitted parameter will be Beta
parg2<-c(Beta=0.5, parg1["Theta"])
# Generate a likelihood map
# [default Phi=seq(from=0.1, to=20, length.out=100) but it is very long]
# Take care, it takes 20 hours ! The data map_Gratiot has the result
\dontrun{
map_Gratiot<-map_phenology(data=data_Gratiot, Phi=seq(from=0.1, to=20, length.out=100),
		parametersfit=parg2, parametersfixed=pfixed)
}
data(map_Gratiot)
# Plot the map
plot(map_Gratiot, col=heat.colors(128))
# Plot the min(-Ln L) for Phi varying at any delta value
plot_phi(map=map_Gratiot)
# Plot the min(-Ln L) for Delta varying with Phi equal to the value for maximum likelihood
plot_delta(map=map_Gratiot)
# Plot the min(-Ln L) for Delta varying with Phi the nearest to 15
plot_delta(map=map_Gratiot, Phi=15)
}
\author{
  Marc Girondot
}

