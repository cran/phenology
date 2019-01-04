#' plot.TaglossData plots formated data used for tagloss analysis
#' @title Plot data used for tagloss analysis.
#' @author Marc Girondot
#' @return Nothing
#' @param x A result for Tagloss_format.
#' @param ... Graphic parameters, see par().
#' @param categories Categories to display.
#' @param col The ramp of colors used for the categories.
#' @param title.legend Title for legend box.
#' @param categories.legend Name of categories to show in legend box.
#' @param show.legend Should the legend box be shown ?
#' @description This function plots the result of Tagloss_format().\cr
#' The default ramp of colors is a grey ramp.
#' @family Model of Tag-loss
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' data_f_21 <- Tagloss_format(outLR, model="21")
#' plot(data_f_21)
#' }
#' @method plot TaglossData
#' @export


plot.TaglossData <- function (x, ..., 
                              categories=c("N22", "N21", "N11", "N10", "N20"), 
                              col=grey(seq(from=0.9, to =0, length.out = length(categories))), 
                              title.legend="Tag history", 
                              categories.legend=categories, 
                              show.legend=TRUE) {
  p3p <- list(...)
  class(x) <- "data.frame"
  par(mar=c(4, 2, 1, 1)+0.4)
  
  p3px <- modifyList(list(xlim=c(0, Tagloss_daymax(x)), ylim=c(1, nrow(x)), 
                  type="n", xlab="Days after first tagging", 
                  ylab="", yaxt="n", bty="n"), p3p)
  p3px <- modifyList(p3px, list(x=1, y=1, type="n", ylab="", col="black"))
  
  do.call(plot, p3px)
  
  i <- 1
  xx <- data.frame(ifelse(is.na(x[, categories[i]]), 0, x[, categories[i]]))
  if (length(categories) > 1) {
    for (i in 2:length(categories)) {
      xx <- cbind(xx, data.frame(ifelse(is.na(x[, categories[i]]), 0, x[, categories[i]])))
    }
  }
  colnames(xx) <- categories
  
  o <- order(rowSums(xx))
  
  x <- xx[o, ]
  
  c <- col
  
  for (r in 1:nrow(x)) {
    
    for (i in 1:ncol(x)) {
      assign(categories[i], value=x[r, categories[i]])
    }
    
    i <- 1
    segments(x0=0, x1=get(categories[i]), y0=r, y1=r, c[i])
    x0 <- get(categories[i])
    
    for (i in 2:length(categories)) {
      if (get(categories[i]) != 0) {
        segments(x0=x0, x1=x0+get(categories[i]), y0=r, y1=r, lty=1, col=c[i])
        x0 <- x0 + get(categories[i])
      }
    }
  }
  
  ylab <- p3p$ylab
  if (is.null(ylab)) ylab <- paste0("Individuals (n=", nrow(x), ")")
  
  mtext(text = ylab, side = 2)
  
  if (show.legend) {
  legend("bottomright", legend = categories.legend, lty=1, 
         col=c, title=title.legend)
  }
  
}
