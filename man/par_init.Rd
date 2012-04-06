\name{par_init}
\alias{par_init}
\title{Calculate initial set of parameters.}
\usage{
  par_init(data = dta, parametersfixed = NULL,
    help = FALSE)
}
\arguments{
  \item{data}{dataset generated with add_format}

  \item{parametersfixed}{Set of fixed parameters}

  \item{help}{If TRUE, an help is displayed}
}
\value{
  The initial set of parameters
}
\description{
  This function is used to generate a first set of
  parameters that is expected to be not to far from the
  final.
}
\details{
  par_init calculates initial set of parameters.
}
\author{
  Marc Girondot
}

