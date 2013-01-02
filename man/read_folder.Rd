\name{read_folder}
\alias{read_folder}
\title{Open all file present in a folder and create a list with these files}
\usage{
  read_folder(folder = ".", read = read.delim, ...)
}
\arguments{
  \item{folder}{Where to search for files}

  \item{read}{Function used to read file. Ex: read.delim}

  \item{...}{Parameters send to the read function}
}
\value{
  Return a list of the data in the files of the folder
  (directory for windows users)
}
\description{
  To create a new list, the syntaxe is \cr
  datalist<-read_folder(folder=".", read=read.delim,
  header=FALSE)\cr
}
\details{
  read_folder opens all files from a folder
}
\examples{
\dontrun{
library(phenology)
# Read all the files from a folder/directory
Gratiot<-read_folder(".")
}
}
\author{
  Marc Girondot
}

