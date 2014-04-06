\name{add_phenology}
\alias{add_phenology}
\title{Create a new dataset or add a timeserie to a previous dataset.}
\usage{
add_phenology(add = file.choose(), name = NULL, reference = NULL,
  month_ref = NULL, header = NULL, format = NULL, previous = NULL,
  silent = FALSE, help = FALSE)
}
\arguments{
  \item{previous}{Name of previous data formated with
  add_phenology or NULL [default] if no previous data
  exists}

  \item{add}{The data to be added. It can be a set of
  several entities that uses the same reference and date
  format}

  \item{name}{The name of the monitored site}

  \item{reference}{as.Date('2001-12-31') The date used as
  1st date}

  \item{month_ref}{If no reference date is given, use this
  month as a reference}

  \item{header}{If the data is read from a file, can be
  used to force header or not}

  \item{format}{The format of the date in the file. Several
  format can be set and the last one that give compatible
  result is used}

  \item{silent}{Does information about added timeseries is
  shown}

  \item{help}{If TRUE, an help is displayed}
}
\value{
Return a list of formated data that can be used ith
fit_phenology()
}
\description{
To create a new dataset, the syntaxe is \cr
data<-add_phenology(add=newdata, name="Site",
reference=as.Date('2001-12-31'), format='%d/%m/%y')\cr To
add a dataset to a previous one, the syntaxe is \cr
data<-add_phenology(previous=previousdata, add=newdata,
name='Site', \cr reference=as.Date('2001-12-31'),
adjust_ref=TRUE, format='%d/%m/%y')\cr To add several
timeseries at the same time with '%d/%m/%y'or '%d/%m/%Y'
date format:\cr data<-add_phenology(add=list(newdata1,
newdata2), name=c('Site1', 'Site2'),\cr
reference=as.Date('2001-12-31'), format=c('%d/%m/%y',
'%d/%m/%Y'))\cr The dataset to be added must include 2 or 3
columns.\cr The first one is the date in the format
specified by\cr the parameter format=. If the number of
nests is known\cr for an exact data, then only one date
must be indicated\cr If the number of nests is known for a
range of date, the\cr first and last dates must be
separated but a - (dash).\cr For example:
1/2/2000-10/2/2000\cr The second column is the number of
nests observed for\cr this date or this range of dates.\cr
The third column is optional and is the name of the
rookery.\cr If only two columns are indicated, the name can
be indicated as\cr a parameter of the function with name=.
If no name is indicated,\cr the default name Site will be
used, but take care, only one \cr rookery of this name can
be used.\cr Several rookeries can be included in the same
file but in this case\cr the rookery name is obligatory at
the third column.\cr The simplest use of this function is
just: \cr phen <- add_phenology()
}
\details{
add_phenology creates a new dataset.
}
\examples{
\dontrun{
library(phenology)
# Read a file with data
Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
data(Gratiot)
# Generate a formatted list nammed data_Gratiot
refdate <- as.Date("2001-01-01")
data_Gratiot<-add_phenology(Gratiot, name="Complete",
	reference=refdate, format="\%d/\%m/\%Y")
# Generate initial points for the optimisation
parg<-par_init(data_Gratiot, parametersfixed=NULL)
# Run the optimisation
result_Gratiot<-fit_phenology(data=data_Gratiot, parametersfit=parg,
	parametersfixed=NULL, trace=1)
data(result_Gratiot)
# Plot the phenology and get some stats
output<-plot(result_Gratiot)
}
}
\author{
Marc Girondot
}

