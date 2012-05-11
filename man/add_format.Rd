\name{add_format}
\alias{add_format}
\title{Create a new dataset or add a timeserie to a previous one.}
\usage{
  add_format(origin = NULL,
    add = stop("A dataset must be indicated"), name = NULL,
    reference = NULL, format = "\%d/\%m/\%y", help = FALSE)
}
\arguments{
  \item{origin}{previousdata or NULL if no previous data
  exists}

  \item{add}{newdata The data to be added. It can be a set
  of several entities that uses the same reference and date
  format}

  \item{name}{'Site' The name of the monitored site}

  \item{reference}{as.Date('2001-12-31') The date used as
  1st date}

  \item{format}{The format of the date in the file. Several
  format can be set and the last one that give compatible
  result is used}

  \item{help}{If TRUE, an help is displayed}
}
\value{
  Return a list of formated data
}
\description{
  To create a new dataset, the syntaxe is \cr
  data<-add_format(add=newdata, name="Site",
  reference=as.Date('2001-12-31'), format='%d/%m/%y')\cr To
  add a dataset to a previous one, the syntaxe is \cr
  data<-add_format(origin=previousdata, add=newdata,
  name='Site', reference=as.Date('2001-12-31'),
  format='%d/%m/%y')\cr To add several timeseries at the
  same time with '%d/%m/%y'or '%d/%m/%Y' date format:
  data<-add_format(add=list(newdata1, newdata2),
  name=c('Site1', 'Site2'),
  reference=as.Date('2001-12-31'), format=c('%d/%m/%y',
  '%d/%m/%Y'))\cr \cr The dataset to be added must include
  2 or 3 columns.\cr The first one is the date in the
  format specified by the parameter format=. If the number
  of nests is known for an exact data, then only one date
  must be indicated If the number of nests is known for a
  range of date, the first and last dates must be separated
  but a - (dash).\cr For example: 1/2/2000-10/2/2000\cr The
  second column is the number of nests observed for this
  date or this range of dates.\cr The third column is
  optional and is the name of the rookery.\cr If only two
  columns are indicated, the name can be indicated as a
  parameter of the function with name=. If no name is
  indicated, the default name Site will be used, but take
  care, only one rookery this name can be used.\cr Several
  rookeries can be included in the same file but in this
  case the rookery name is obligatory at the third column.
}
\details{
  add_format creates a new dataset.
}
\examples{
library(phenology)
# Read a file with data
# Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", , header=FALSE)
data(Gratiot)
# Generate a formatted list nammed data_Gratiot
data_Gratiot<-add_format(origin=NULL, add=Gratiot, name="Complete", reference=as.Date("2001-01-01"), format="\%d/\%m/\%Y")
# Generate initial points for the optimisation
parg<-par_init(data_Gratiot, parametersfixed=NULL)
# Run the optimisation
# result_Gratiot<-fit_phenology(data=data_Gratiot, parametersfit=parg, parametersfixed=NULL, trace=1)
data(result_Gratiot)
# Plot the phenology and get some stats
output<-plot_phenology(result=result_Gratiot, pdf=FALSE)
}
\author{
  Marc Girondot
}

