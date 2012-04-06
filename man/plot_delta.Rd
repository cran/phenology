\name{plot_delta}
\alias{plot_delta}
\title{Plot a likelihood lineplot obtained after map_phenology.}
\usage{
  plot_delta(map = NULL, Phi = NULL, pdf = FALSE,
    pdfname = "Map.pdf", help = FALSE)
}
\arguments{
  \item{map}{A map generated with map_phenology}

  \item{Phi}{Phi value or NULL}

  \item{pdf}{TRUE or FALSE, indicates if a pdf file is
  generated.}

  \item{pdfname}{Name of pdf file}

  \item{help}{If TRUE, an help is displayed}
}
\value{
  Return None
}
\description{
  This function plots a likelihood lineplot obtained after
  map_phenology.
}
\details{
  plot_delta plots the likelihood delta for fixed Phi
  value.
}
\author{
  Marc Girondot
}

