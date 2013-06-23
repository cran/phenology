#' ScalePreviousPlot returns the scale of the previous plot
#' @title Return the scale of the previous plot
#' @author Marc Girondot
#' @return A list with xlim and ylim
#' @description Return a list with the limits of the previous plot
#' @examples
#' par(xaxs="i", yaxs="i")
#' plot(x=1:100, y=sin(1:100), type="l", bty="n", xlim=c(1,200), xlab="x", ylab="y")
#' xlim= ScalePreviousPlot()$xlim
#' ylim= ScalePreviousPlot()$ylim
#' par(xaxs="r", yaxs="i")
#' plot(x=1:100, y=sin(1:100), type="l", bty="n", xlim=c(1,200), xlab="x", ylab="y")
#' xlim= ScalePreviousPlot()$xlim
#' ylim= ScalePreviousPlot()$ylim
#' @export


ScalePreviousPlot <- function() {
  if (par("xaxs")=="i") {
    x1 <- par("usr")[1]
    x2 <- par("usr")[2]
  } else {
    x2 <- (par("usr")[1]+par("usr")[2]*26)/27
    x1 <- x2*26-par("usr")[2]/0.04
  }
    if (par("yaxs")=="i") {
    y1 <- par("usr")[3]
    y2 <- par("usr")[4]
  } else {
    y2 <- (par("usr")[3]+par("usr")[4]*26)/27
    y1 <- y2*26-par("usr")[4]/0.04
  }
  return(list(xlim=c(x1, x2), ylim=c(y1, y2)))
}