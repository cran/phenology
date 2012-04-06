\name{map_phenology}
\alias{map_phenology}
\title{Generate a likelihood map varying Phi and Delta.}
\usage{
  map_phenology(data = NULL, parametersfit = NULL,
    parametersfixed = NA,
    Phi = seq(from = 0.2, to = 20, length.out = 100),
    Delta = NULL, method_incertitude = 2,
    zero_counts = TRUE, help = FALSE)
}
\arguments{
  \item{data}{dataset generated with add_format}

  \item{parametersfixed}{Set of fixed parameters}

  \item{parametersfit}{Set of parameters to be fitted}

  \item{Phi}{Phi values to be analyzed}

  \item{Delta}{Delta value to be analyzed}

  \item{method_incertitude}{2 [default] is the correct one
  from a statistical point of view; 0 is an aproximate
  method more rapid; 1 is an alternative more rapid but
  biased.}

  \item{zero_counts}{Example c(TRUE, TRUE, FALSE) indicates
  whether the zeros have been recorder for each of these
  timeseries. Defaut is TRUE for all.}

  \item{help}{If TRUE, an help is displayed}
}
\value{
  Display a likelihood map
}
\description{
  This function generates a map of likelihood varying Phi
  and Delta.
}
\details{
  map_phenology generates a likelihood map.
}
\author{
  Marc Girondot
}

