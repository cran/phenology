\name{plot_map}
\alias{plot_map}
\title{Plot a likelihood map with Delta and Phi varying.}
\usage{
  plot_map(map = NULL, pdf = FALSE, pdfname = "Map.pdf",
    col = heat.colors(128), help = FALSE)
}
\arguments{
  \item{map}{A map generated with map_phenology.}

  \item{pdf}{TRUE or FALSE, indicates if a pdf file is
  generated.}

  \item{pdfname}{Name of file if pdf=TRUE}

  \item{col}{Colors could be heat.colors(128) or
  rainbow(64) or col=gray(c(seq(0, 1, length.out=128)))}

  \item{help}{If TRUE, an help is displayed}
}
\value{
  Return None
}
\description{
  This function plots a likelihood map obtained after
  map_phenology.
}
\details{
  plot_map plots a likelihood map with Delta and Phi
  varying.
}
\author{
  Marc Girondot
}

