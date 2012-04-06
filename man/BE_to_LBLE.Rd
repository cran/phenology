\name{BE_to_LBLE}
\alias{BE_to_LBLE}
\title{Transform a set of parameters from Begin End format to LengthB LengthE.}
\usage{
  BE_to_LBLE(parameters = NULL, help = FALSE)
}
\arguments{
  \item{parameters}{Set of current parameters}

  \item{help}{If TRUE, an help is displayed}
}
\value{
  Return the set of modified parameters
}
\description{
  This function is used to transform a set of parameters
  that uses Begin, Peak and End to a set of parameters that
  uses LengthB, Peak and LengthE.
}
\details{
  BE_to_LBLE transforms a set of parameters from Begin End
  format to LengthB LengthE.
}
\author{
  Marc Girondot
}

