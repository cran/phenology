\name{plot_add}
\alias{plot_add}
\title{Add a plot to a previous one}
\usage{
plot_add(...)
}
\arguments{
  \item{...}{Parameters for plot()}
}
\value{
Nothing
}
\description{
To plot data, just add use it as a normal plot. It will
plot the new data without axes, or labels for axes.
}
\details{
plot_add adds a plot to a previous one
}
\examples{
plot(x=1:100, y=sin(1:100), type="l", bty="n", xlim=c(1,200), xlab="x", ylab="y")
plot_add(x=1:200, y=cos(1:200), type="l", bty="n", col="red")
}
\author{
Marc Girondot
}

