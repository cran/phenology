#' plot.Remigration plots the remigration intervals.
#' @title Plot the remigration intervals.
#' @author Marc Girondot
#' @return An invisible dataframe with values used for plotting.
#' @param x Object obtained from Bayesian.remigration()
#' @param legend TRUE or FALSE or c(x, y)
#' @param ... Parameters transmitted to plot
#' @description Plot the remigration intervals.\cr
#' @family Model of Remigration Interval
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' 
#' # Each year a fraction of 0.9 is surviving
#' s <- c(s=0.9)
#' # Probability of tag retention; 0.8
#' t <- c(t=0.8)
#' # Time-conditional return probability - This is the true remigration rate
#' r <- c(r1=0.1, r2=0.8, r3=0.7, r4=0.7, r5=1)
#' # Capture probability
#' p <- c(p1=0.6, p2=0.6, p3=0.6, p4=0.6, p5=0.6)
#' # Number of observations for 400 tagged females after 1, 2, 3, 4, and 5 years
#' OBS <- c(400, 10, 120, 40, 20, 10)
#' 
#' kl_s <- length(s)
#' kl_t <- length(t)
#' kl_r <- length(r)
#' kl_p <- length(p)
#' 
#' pMCMC <- data.frame(Density=c("newdbeta", "newdbeta", rep("dunif", kl_r), 
#'                               rep("newdbeta", kl_p), "dunif"), 
#'                     Prior1=c(s, t, rep(0, kl_r), rep(0.2, kl_p), 0), 
#'                     Prior2=c(0.02, 0.02, rep(1, kl_r), rep(0.08, kl_p), 10), 
#'                     SDProp=c(0.05, 0.05, rep(0.05, kl_r), rep(0.05, kl_p), 0.05), 
#'                  Min=c(0, 0, rep(0, kl_r), rep(0, kl_p), 0),  
#'                  Max=c(1, 1, rep(1, kl_r), rep(1, kl_p), 10),  
#'                  Init=c(s, t, r, p, 1), stringsAsFactors = FALSE, 
#'                  row.names=c("s", 
#'                                 "t", 
#'                                 names(r), 
#'                                 names(p), "sd")
#' )
#' rMCMC <- Bayesian.remigration(parameters = pMCMC, 
#' n.iter = 1000000, 
#' n.adapt = 300000,
#' trace=10000, 
#' data=OBS)
#' 
#' plot(rMCMC)
#' 
#' }
#' @method plot Remigration
#' @export


plot.Remigration <- function(x, legend=TRUE, ...) {
  ri <- x$RI
  p3p <- list(...)
  p3p <- modifyList(list(x=ri, freq=FALSE, las=1, main="", 
                         xlab="Remigration interval"), 
                    p3p)
  do.call(hist, p3p)
  if (legend != FALSE) {
    if (is.numeric(legend)) {
      posx <- legend
    } else {
      posx <- ScalePreviousPlot()$xlim["begin"]+ScalePreviousPlot()$xlim["range"]*3/5
    }
    text(x=posx, y=ScalePreviousPlot()$ylim["end"]*2/3, labels=paste0("Mean: ", format(mean(ri), digits=3)), pos=4)
    text(x=posx, y=ScalePreviousPlot()$ylim["end"]*1.7/3, labels=paste0("SE: ", format(sd(ri), digits=2)), pos=4)
    text(x=posx, y=ScalePreviousPlot()$ylim["end"]*1.4/3, labels=paste0("Quantile 2.5%: ", format(quantile(ri, probs = 0.025), digits=3)), pos=4)
    text(x=posx, y=ScalePreviousPlot()$ylim["end"]*1.1/3, labels=paste0("Quantile 97.5%: ", format(quantile(ri, probs = 0.975), digits=3)), pos=4)
  }
  par(xpd=TRUE)
  segments(x0=mean(ri), y0=0, y1=ScalePreviousPlot()$ylim["end"]*1.1, lty=4)
  
  segments(x0=quantile(ri, probs = 0.025), y0=ScalePreviousPlot()$ylim["end"]*1.1, 
           x1=quantile(ri, probs = 0.975), lty=1)
  segments(x0=quantile(ri, probs = 0.025), y0=ScalePreviousPlot()$ylim["end"]*1.05, 
           y1=ScalePreviousPlot()$ylim["end"]*1.15, lty=1)
  segments(x0=quantile(ri, probs = 0.975), y0=ScalePreviousPlot()$ylim["end"]*1.05, 
           y1=ScalePreviousPlot()$ylim["end"]*1.15, lty=1)
}
