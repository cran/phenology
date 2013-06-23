\name{plot_errbar}
\alias{plot_errbar}
\title{Plot a xy graph with error bar on x and/or y}
\usage{
  plot_errbar(..., errbar.x = NULL, errbar.y = NULL,
    errbar.x.plus = NULL, errbar.x.minus = NULL,
    errbar.y.plus = NULL, errbar.y.minus = NULL,
    errbar.tick = 1/50, errbar.lwd = par("lwd"),
    errbar.lty = par("lty"), errbar.col = par("fg"))
}
\arguments{
  \item{...}{Parameters for plot() such as main= or ylim=}

  \item{errbar.x}{The length of error bars for x. Recycled
  if necessary.}

  \item{errbar.x.plus}{The length of positive error bars
  for x. Recycled if necessary.}

  \item{errbar.x.minus}{The length of negative error bars
  for x. Recycled if necessary.}

  \item{errbar.y}{The length of error bars for y. Recycled
  if necessary.}

  \item{errbar.y.plus}{The length of positive error bars
  for y. Recycled if necessary.}

  \item{errbar.y.minus}{The length of negative error bars
  for y. Recycled if necessary.}

  \item{errbar.tick}{Size of small ticks at the end of
  error bars defined as a proportion of total width or
  height graph size.}

  \item{errbar.lwd}{Error bar line width, see par("lwd")}

  \item{errbar.lty}{Error bar line type, see par("lwd")}

  \item{errbar.col}{Error bar line color, see par("col")}
}
\value{
  Nothing
}
\description{
  To plot data, just add use it as a normal plot but add
  the errbar.x and errbar.y values.
}
\details{
  plot_errbar plot a xy graph with error bar on x and/or y
}
\examples{
plot_errbar(1:100, rnorm(100, 1, 2),
		xlab="axe x", ylab="axe y", bty="n", xlim=c(1,100),
		errbar.x=2, errbar.y=rnorm(100, 1, 0.1))
}
\author{
  Marc Girondot
}

