\name{LBLE_to_BE}
\alias{LBLE_to_BE}
\title{Transform a set of parameters from LengthB LengthE to Begin End format.}
\usage{
  LBLE_to_BE(parameters = NULL, help = FALSE)
}
\arguments{
  \item{parameters}{Set of current parameters}

  \item{help}{If TRUE, an help is displayed}
}
\value{
  Return a set of modified parameters
}
\description{
  This function is used to transform a set of parameters
  that uses LengthB, Peak and LengthE to a set of
  parameters that uses Begin, Peak and End.
}
\details{
  LBLE_to_BE transforms a set of parameters from LengthB
  LengthE to Begin End format.
}
\author{
  Marc Girondot
}

