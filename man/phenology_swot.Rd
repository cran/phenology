\name{phenology_swot}
\alias{phenology_swot}
\title{Fit a nesting season of marine turtles}
\usage{
phenology_swot(header = NULL, reference = NULL, month_ref = NULL,
  format = NULL)
}
\arguments{
  \item{header}{Does the timeseries has header.}

  \item{reference}{Date used as reference. Is the day 1.}

  \item{month_ref}{Reference month. Generally will be 1
  (January) or 7 (July).}

  \item{format}{Format for dates.}
}
\value{
Return a list of with data and result
}
\description{
This function tries to make very simple the use of this
package.\cr In most of the cases, values for header,
reference and format need not to be set\cr because they are
detected automatically.
}
\details{
phenology_swot is a simplified function for phenology
package developped as part of SWOT project
}
\author{
Marc Girondot \email{marc.girondot@u-psud.fr}
}

