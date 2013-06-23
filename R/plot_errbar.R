#' plot_errbar plot a xy graph with error bar on x and/or y
#' @title Plot a xy graph with error bar on x and/or y
#' @author Marc Girondot
#' @return Nothing
#' @param ... Parameters for plot() such as main= or ylim=
#' @param errbar.x The length of error bars for x. Recycled if necessary.
#' @param errbar.x.plus The length of positive error bars for x. Recycled if necessary.
#' @param errbar.x.minus The length of negative error bars for x. Recycled if necessary.
#' @param errbar.y The length of error bars for y. Recycled if necessary.
#' @param errbar.y.plus The length of positive error bars for y. Recycled if necessary.
#' @param errbar.y.minus The length of negative error bars for y. Recycled if necessary.
#' @param errbar.tick Size of small ticks at the end of error bars defined as a proportion of total width or height graph size.
#' @param errbar.lwd Error bar line width, see par("lwd")
#' @param errbar.lty Error bar line type, see par("lwd")
#' @param errbar.col Error bar line color, see par("col")
#' @description To plot data, just add use it as a normal plot but add the errbar.x and errbar.y values.
#' @examples
#' plot_errbar(1:100, rnorm(100, 1, 2), 
#'		xlab="axe x", ylab="axe y", bty="n", xlim=c(1,100), 
#' 		errbar.x=2, errbar.y=rnorm(100, 1, 0.1))
#' @export


plot_errbar <- function(..., 
                        errbar.x=NULL, errbar.y=NULL, 
                        errbar.x.plus=NULL, errbar.x.minus=NULL, 
                        errbar.y.plus=NULL, errbar.y.minus=NULL, 
                        errbar.tick=1/50, 
                        errbar.lwd=par("lwd"), 
                        errbar.lty=par("lty"), 
                        errbar.col=par("fg")) 
  {
  par.plot <- list(...)
  do.call(plot, par.plot) 
  
  x <- par.plot[["x"]]
  if (is.null(x)) x <- par.plot[[1]]
  if (is.data.frame(x) | is.matrix(x)) {
    y <- x[,2]
    x <- x[,1]
  } else {
    y <- par.plot[["y"]]
  }
  if (is.null(y)) y <- par.plot[[2]]
  
  if (is.null(errbar.x.minus) & !is.null(errbar.x)) {
  	errbar.x.minus <- errbar.x
  }
  if (is.null(errbar.x.plus) & !is.null(errbar.x)) {
  	errbar.x.plus <- errbar.x
  }
  if (is.null(errbar.y.minus) & !is.null(errbar.y)) {
  	errbar.y.minus <- errbar.y
  }
  if (is.null(errbar.y.plus) & !is.null(errbar.y)) {
  	errbar.y.plus <- errbar.y
  }

  sizebar <- (par("usr")[4]-par("usr")[3])*errbar.tick

  if (!is.null(errbar.x.minus)) {
    segments(x-errbar.x.minus, y, x, y, 
             col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
    segments(x-errbar.x.minus, y-sizebar, x-errbar.x.minus, y+sizebar, 
             col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
  }
  if (!is.null(errbar.x.plus)) {
    segments(x+errbar.x.plus, y, x, y, 
             col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
    segments(x+errbar.x.plus, y-sizebar, x+errbar.x.plus, y+sizebar, 
             col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
  }
    
  sizebar <- (par("usr")[2]-par("usr")[1])*errbar.tick
  
  if (!is.null(errbar.y.minus)) {
    segments(x, y-errbar.y.minus, x, y, 
             col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
    segments(x-sizebar, y-errbar.y.minus, x+sizebar, y-errbar.y.minus, 
             col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
  }
  if (!is.null(errbar.y.plus)) {
    segments(x, y+errbar.y.plus, x, y, 
             col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
    segments(x-sizebar, y+errbar.y.plus, x+sizebar, y+errbar.y.plus, 
             col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
  }

}
