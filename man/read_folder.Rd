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
library(phenology)
# Read all the files from a folder/directory
# Gratiot<-read_folder(".")
data(Gratiot)
# Generate a formatted list nammed data_Gratiot
data_Gratiot<-add_format(add=Gratiot, name="Complete", reference=as.Date("2001-01-01"), format="\%d/\%m/\%Y")
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

