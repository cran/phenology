\name{MinBMinE_to_Min}
\alias{MinBMinE_to_Min}
\title{Transform a set of parameters from MinB and MinE to Min}
\usage{
MinBMinE_to_Min(parameters = stop("A set of parameters must be indicated"))
}
\arguments{
  \item{parameters}{Set of current parameters}
}
\value{
Return a set of modified parameters
}
\description{
This function is used to transform a set of parameters that
uses MinB and MinE to a set of parameters that uses Min.
}
\details{
MinBMinE_to_Min transforms a set of parameters from MinB
and MinE to Min
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
parg1<-MinBMinE_to_Min(parameters=parg)
}
\author{
Marc Girondot
}

