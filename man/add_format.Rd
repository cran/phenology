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

  \item{add}{newdata The data to be added}

  \item{name}{'Site' The name of the monitored site}

  \item{reference}{as.Date('2001-12-31') The date used as
  1st date}

  \item{format}{The format of the date in the file}

  \item{help}{If TRUE, an help is displayed}
}
\value{
  Return a list of formated data
}
\description{
  To create a new dataset, the syntaxe is
  data<-add_format(add=newdata, name="Site", +
  reference=as.Date('2001-12-31'), format='%d/%m/%y') To
  add a dataset to a previous one, the syntaxe is
  data<-add_format(origin=previousdata, add=newdata,
  name='Site', + reference=as.Date('2001-12-31'),
  format='%d/%m/%y')

  The dataset to be added must include 2 or 3 columns. The
  first one is the date in the format specified by the
  parameter format=. If the number of nests is known for an
  exact data, then only one date must be indicated If the
  number of nests is known for a range of date, the first
  and last dates must be separated but a -. For example:
  1/2/2000-10/2/2000 The second column is the number of
  nests observed for this date or this range of dates. The
  third column is optional and is the name of the rookery.
  If only two columns are indicated, the name can be
  indicated as a parameter of the function with name=. If
  no name is indicated, the default name Site will be used,
  but take care, only one rookery this name can be used.
  Several rookeries can be included in the same file but in
  this case the rookery name is obligatory at the third
  column.
}
\details{
  add_format creates a new dataset.
}
\author{
  Marc Girondot
}

