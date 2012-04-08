\name{plot_phi}
\alias{plot_phi}
\title{Plot the best likelihood for fixed Phi value.}
\usage{
  plot_phi(map = NULL, pdf = FALSE, pdfname = "Map.pdf",
    help = FALSE)
}
\arguments{
  \item{map}{A map generated with map_phenology}

  \item{pdf}{TRUE or FALSE, indicates if a pdf file is
  generated.}

  \item{pdfname}{Name of pdf file}

  \item{help}{If TRUE, an help is displayed}
}
\value{
  Return None
}
\description{
  The function "plot_phi" plots the best likelihood for
  each Phi value.
}
\details{
  plot_phi plots the best likelihood for fixed Phi value.
}
\examples{
library("phenology")
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
# Extract the fitted parameters
parg1<-extract_result(result_Gratiot)
# Add constant Alpha and Tau values [day d amplitude=(Alpha+Nd*Beta)^Tau with Nd being the number of counts for day d]
pfixed<-c(parg1, Alpha=0, Tau=1)
pfixed<-pfixed[-which(names(pfixed)=="Theta")]
# The only fitted parameter will be Beta
parg2<-c(Beta=0.5, parg1["Theta"])
# Generate a likelihood map [default Phi=seq(from=0.1, to=20, length.out=100) but it is very long]
# Take care, it takes 20 hours ! The data map_Gratiot has the result
# map_Gratiot<-map_phenology(data=data_Gratiot, Phi=seq(from=0.1, to=20, length.out=100), parametersfit=parg2, parametersfixed=pfixed)
data(map_Gratiot)
# Plot the min(-Ln L) for Phi varying at any delta value
plot_phi(map=map_Gratiot, pdf=FALSE)
}
\author{
  Marc Girondot
}

