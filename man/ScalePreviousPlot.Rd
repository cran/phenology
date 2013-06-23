\name{ScalePreviousPlot}
\alias{ScalePreviousPlot}
\title{Return the scale of the previous plot}
\usage{
  ScalePreviousPlot()
}
\value{
  A list with xlim and ylim
}
\description{
  Return a list with the limits of the previous plot
}
\details{
  ScalePreviousPlot returns the scale of the previous plot
}
\examples{
par(xaxs="i", yaxs="i")
plot(x=1:100, y=sin(1:100), type="l", bty="n", xlim=c(1,200), xlab="x", ylab="y")
xlim= ScalePreviousPlot()$xlim
ylim= ScalePreviousPlot()$ylim
par(xaxs="r", yaxs="i")
plot(x=1:100, y=sin(1:100), type="l", bty="n", xlim=c(1,200), xlab="x", ylab="y")
xlim= ScalePreviousPlot()$xlim
ylim= ScalePreviousPlot()$ylim
}
\author{
  Marc Girondot
}

