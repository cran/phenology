#' plot.TableECFOCF plots a TableECFOCF dataset.
#' @title Plot a TableECFOCF dataset.
#' @author Marc Girondot
#' @return Nothing
#' @param x A CMR file summarized using TableECFOCF()
#' @param ... Graphic parameters
#' @param result What should be plotted: ECFOCF or data, ECF, OCF.
#' @param period The period that will be plotted.
#' @param cex.points The maximum magnification to be used for points relative to the current setting of cex.
#' @param pch Character to be used for points.
#' @param col Color to be used for points.
#' @param cex.axis The magnification to be used for axis annotation relative to the current setting of cex.
#' @param show.labels Logical to be used to show figures.
#' @param col.labels Color of figures.
#' @param cex.labels The magnification to be used for figures.
#' @param show.0 Logical to show 0 counts.
#' @param pch.0 Character used for 0 counts.
#' @param cex.0 The magnification to be used for character for 0 counts.
#' @param col.0 Color of characters for 0 counts.
#' @param show.scale If TRUE, show the scale as a legend
#' @description This function plots a CMR file summarized using TableECFOCF().\cr
#' @family Model of Clutch Frequency
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' data(MarineTurtles_2002)
#' ECFOCF_2002 <- TableECFOCF(MarineTurtles_2002)
#' par(mar=c(4, 4, 1, 1)+0.4)
#' plot(ECFOCF_2002, bty="n", las=1, cex.points=3, 
#'      cex.axis = 0.8, main="Year 2002")
#' plot(ECFOCF_2002, bty="n", las=1, cex.points=5, cex.0=0.2, 
#'      col="red", show.0 = TRUE, col.0="blue")
#' plot(ECFOCF_2002, bty="n", las=1, cex.points=3, col="lightgrey",  
#'      col.labels = "red", show.labels=TRUE)
#' plot(ECFOCF_2002, bty="n", las=1, cex.points=3, pch=NA, 
#'      col.labels = "red", show.labels=TRUE)
#' plot(ECFOCF_2002, bty="n", las=1, cex.points=3, pch=NA, 
#'      col.labels = "red", show.labels=TRUE, cex.0=0.2, 
#'      show.0 = TRUE, col.0="blue", pch.0=4)
#' plot(ECFOCF_2002, bty="n", las=1, result="OCF")
#' plot(ECFOCF_2002, bty="n", las=1, result="ECF")
#' plot(ECFOCF_2002, bty="n", las=1, result="ECF", type="l", main="2002 season", 
#'      xlab="Clutch frequency")
#' par(new=TRUE)
#' plot(ECFOCF_2002, bty="n", las=1, result="OCF", type="l", main="", 
#'      ylim=ScalePreviousPlot()$ylim[c("begin", "end")], 
#'      xlab="", ylab="", 
#'      col="red", 
#'      xaxt="n", yaxt="n", axes=FALSE)
#' legend("topright", legend=c("OCF", "ECF"), lty=1, col=c("red", "black"))
#' 
#' ECFOCF_2002 <- TableECFOCF(MarineTurtles_2002, date0=as.Date("2002-01-01"))
#' 
#' plot(ECFOCF_2002, period=13)
#' }
#' @method plot TableECFOCF
#' @export


# plot de la table ECF OCF ####

plot.TableECFOCF <- function(x, ..., result="ecfocf", 
                             period=1, 
                             cex.points=4, 
                             pch=19, 
                             col="black",
                             cex.axis=0.8, 
                             cex.labels=0.5, 
                             col.labels="red", 
                             show.labels=FALSE, 
                             show.0=FALSE, 
                             pch.0=4, 
                             cex.0=0.5, 
                             col.0="blue", 
                             show.scale = TRUE) {
  
  # result="ecfocf";
  # period=1;
  # cex.points=4; 
  # pch=19;
  # col="black";
  # cex.axis=0.8; 
  # cex.labels=0.5; 
  # col.labels="red"; 
  # show.labels=FALSE; 
  # show.0=FALSE; 
  # pch.0=4; 
  # cex.0=0.5;
  # col.0="blue"; 
  # show.scale = TRUE
  
  p3p <- list(...)
  result <- tolower(result)
  
  x <- x[, , period]
  
  if ((result=="ecfocf") | (result == "data")) {
    do.call(plot, modifyList(list(x=1, y=1, type="n", xlim=c(0, ncol(x)-1), 
                                  ylim=c(0, nrow(x)-1), 
                                  xlab="Observed Clutch Frequency", 
                                  ylab="Estimated Clutch Frequency", xaxt="n", yaxt="n"), p3p))
    axis(side=1, at=0:(ncol(x)-1), cex.axis=cex.axis)
    do.call(axis, modifyList(list(side=2, at=0:(nrow(x)-1), cex.axis=cex.axis), p3p)[c("side", "at", "cex.axis", "las", "labels")])
    sc <- max(x, na.rm= TRUE)
    
    for (c in 1:ncol(x)) {
      for (r in 1:nrow(x)) {
        if ((x[r, c] != 0)) {
          do.call(points, 
                  modifyList(list(x=r-1, y=c-1, cex=x[r, c]/sc*cex.points, pch=pch, col=col), 
                             p3p)[c("x", "y", "cex", "pch", "col")])
        }
      }
    }
    
    for (c in 1:ncol(x)) {
      for (r in 1:nrow(x)) {
        if ((x[r, c] == 0) & (show.0)) {
          do.call(points, 
                  modifyList(list(x=r-1, y=c-1, cex=cex.0, pch=pch.0, col=col.0), 
                             p3p)[c("x", "y", "cex", "pch", "col")])
        }
        if (show.labels & (x[r, c] != 0)) {
          do.call(text, 
                  modifyList(list(x=r-1, y=c-1, cex=cex.labels, col=col.labels, labels=as.character(x[r, c])), 
                             p3p)[c("x", "y", "cex", "col", "labels")])
        }
      }
    }
    if (show.scale) {
      legend(x = "bottomright", 
             pch=c(ifelse(show.0, pch.0, NA), pch, pch, pch), 
             col=c(col.0, col, col, col), 
             pt.cex=c(cex.0, (seq(from=0, to=sc, length.out = 4)[2:4])/sc*cex.points), 
             legend=c(0, floor(seq(from=0, to=sc, length.out = 4)[2:4])), 
             title="Scale"
      )
      
      
    }
    
  } else {
    if (result == "ocf") {
      ecfocf <- x
      main="Observed OCF"
      ocf <- rowSums(ecfocf, na.rm = TRUE)
      do.call(plot, modifyList(list(x=0:(length(ocf)-1), 
                                    xlab="Observed Clutch Frequency", 
                                    ylab="Frequency", 
                                    main=main,
                                    col=col, 
                                    cex.axis=cex.axis, 
                                    pch=pch, 
                                    y=ocf, type="h"), p3p)[c("x", "y", "type", "col", 
                                                             "main", "cex.axis", "bty", "las", 
                                                             "xlab", "ylab", "xlim", "ylim", 
                                                             "xaxt", "yaxt", "axes")])
      
    }
    
    if (result == "ecf") {
      ecfocf <- x
      main="Observed ECF"
      ecf <- colSums(ecfocf, na.rm = TRUE)
      do.call(plot, modifyList(list(x=0:(length(ecf)-1), 
                                    xlab="Estimated Clutch Frequency", 
                                    ylab="Frequency", 
                                    main=main, 
                                    col=col,
                                    cex.axis=cex.axis, 
                                    pch=pch, 
                                    y=ecf, type="h"), p3p)[c("x", "y", "type", "col", 
                                                             "main", "cex.axis", "bty", "las", 
                                                             "xlab", "ylab", "xlim", "ylim", 
                                                             "xaxt", "yaxt", "axes")])
      
      
    }
    
    
    
  }
  
  
}

