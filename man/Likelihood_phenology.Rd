\name{Likelihood_phenology}
\alias{Likelihood_phenology}
\title{Estimate the likelihood of timeseries based on a set of parameters.}
\usage{
  Likelihood_phenology(data = dta, parametersfit = x,
    parametersfixed = NULL, method_incertitude = 0,
    help = FALSE)
}
\arguments{
  \item{data}{=dataset generated with add_format}

  \item{parametersfixed}{Set of fixed parameters}

  \item{parametersfit}{Set of parameters to be fitted}

  \item{method_incertitude}{2 [default] is the correct one
  from a statistical point of view; 0 is an aproximate
  method more rapid; 1 is an alternative more rapid but
  biased.}

  \item{help}{If TRUE, an help is displayed}
}
\value{
  Return a list of formated data
}
\description{
  This function is used to estimate the likelihood based on
  a set of parameters.
}
\details{
  Likelihood_phenology estimate likelihood for a set of
  parameters.
}
\author{
  Marc Girondot
}

