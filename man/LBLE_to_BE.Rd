\name{LBLE_to_BE}
\alias{LBLE_to_BE}
\title{Transform a set of parameters from LengthB LengthE to Begin End.}
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
  LengthE to Begin End.
}
\examples{
# Read a file with data
# Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", , header=FALSE)
data(Gratiot)
# Generate a formatted list nammed data_Gratiot
data_Gratiot<-add_format(origin=NULL, add=Gratiot, name="Complete", reference=as.Date("2001-01-01"), format="\%d/\%m/\%Y")
# Generate initial points for the optimisation
parg<-par_init(data_Gratiot, parametersfixed=NULL)
# Change the parameters to Begin End format
parg1<-LBLE_to_BE(parameters=parg)
# And change back to LengthB LengthE
parg2<-BE_to_LBLE(parameters=parg1)
}
\author{
  Marc Girondot
}

