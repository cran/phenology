\name{toggle_Min_PMin}
\alias{toggle_Min_PMin}
\title{Transform a set of parameters from Min, MinB or MinE to PMin, PminB or PminE, or reverse}
\usage{
  toggle_Min_PMin(parameters = stop("A set of parameters must be indicated"))
}
\arguments{
  \item{parameters}{Set of current parameters}
}
\value{
  Return a set of modified parameters
}
\description{
  This function is used to transform a set of parameters
  that uses Min, MinB or MinE to a set of parameters that
  uses PMin, PminB or PminE, or reverse.
}
\details{
  toggle_Min_PMin transforms a set of parameters from Min,
  MinB or MinE to PMin, PminB or PminE, or reverse
}
\examples{
# Read a file with data
# Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", , header=FALSE)
data(Gratiot)
# Generate a formatted list nammed data_Gratiot
refdate <- as.Date("2001-01-01")
data_Gratiot<-add_phenology(Gratiot, name="Complete", reference=refdate, format="\%d/\%m/\%Y")
# Generate initial points for the optimisation
parg<-par_init(data_Gratiot)
# Change the parameters to PMinB and PMinE
parg1<-toggle_Min_PMin(parameters=parg)
# And change back to MinB and MinE
parg2<-toggle_Min_PMin(parameters=parg1)
}
\author{
  Marc Girondot
}

