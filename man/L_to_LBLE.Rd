\name{L_to_LBLE}
\alias{L_to_LBLE}
\title{Transform a set of parameters from Length format to LengthB LengthE}
\usage{
  L_to_LBLE(parameters = stop("Set of parameters must be given"))
}
\arguments{
  \item{parameters}{Set of current parameters}
}
\value{
  Return the set of modified parameters
}
\description{
  This function is used to transform a set of parameters
  that uses Length to a set of parameters uses LengthB and
  LengthE.
}
\details{
  L_to_LBLE transforms a set of parameters from Length
  format to LengthB LengthE.
}
\examples{
\dontrun{
# Read a file with data
Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
}
data(Gratiot)
# Generate a formatted list nammed data_Gratiot
refdate <- as.Date("2001-01-01")
data_Gratiot<-add_phenology(Gratiot, name="Complete",
		reference=refdate, format="\%d/\%m/\%Y")
# Generate initial points for the optimisation
parg<-par_init(data_Gratiot, parametersfixed=NULL)
# Change the parameters to Begin End format
parg1<-LBLE_to_L(parameters=parg)
# And change back to LengthB LengthE.
parg2<-L_to_LBLE(parameters=parg1)
}
\author{
  Marc Girondot
}

