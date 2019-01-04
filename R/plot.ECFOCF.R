#' plot.ECFOCF plots a result of clutch frequency fit.
#' @title Plot a result of clutch frequency fit.
#' @author Marc Girondot
#' @return Nothing
#' @param x A result for fitCF().
#' @param ... Graphic parameters, see plot.TableECFOCF() or par.
#' @param result What result will be plotted: data, dataOCF, dataECF, ECF, OCF, ECFOCF, ECFOCF0, CF, Prob, period
#' @param category What category will be plotted, numeric or NA for all.
#' @param period The period that will be plotted.
#' @description This function plots the result of fitCF().\cr
#' The result \code{data} plots the observed ECF-OCF table.\cr
#' The result \code{dataOCF} plots the observed OCF table.\cr
#' The result \code{dataECF} plots the observed ECF table.\cr
#' The result \code{CF} plots the true clutch frequency.\cr
#' The result \code{OCF} plots the observed clutch frequency.\cr
#' The result \code{ECF} plots the estimated clutch frequency.\cr
#' The result \code{ECFOCF} plots the bivariate observed vs. estimated clutch frequency.\cr
#' The result \code{ECFOCF0} plots the bivariate observed vs. estimated clutch frequency without the 0 OCF.\cr
#' The result \code{prob} plots the probabilities of capture.\cr
#' The result \code{period} plots the probabilities of nesting according to period.\cr
#' If category is left to NA, the compound value for all the population is plotted.\cr
#' When result="data" is used, this is a parser for plot.TableECFOCF().\cr
#' See this function for the parameters.\cr
#' The parameter y.axis is the shift of the x legends for result="prob".
#' @family Model of Clutch Frequency
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' data(MarineTurtles_2002)
#' ECFOCF_2002 <- TableECFOCF(MarineTurtles_2002)
#' o_mu1p2_NB <- fitCF(x = c(mu = 4.6426989650675701, 
#'                          sd = 75.828239144717074, 
#'                          p1 = 0.62036295627161053,
#'                          p2 = -2.3923021862881511, 
#'                          OTN = 0.33107456308054345),
#'                  data=ECFOCF_2002)
#'                  
#' par(mar=c(4, 4, 1, 1)+0.4)
#' plot(o_mu1p2_NB, result="data", category=NA, 
#'      bty="n", las=1, cex.points=3, cex.axis = 0.8)
#' plot(o_mu1p2_NB,result="data", category=NA, 
#'      bty="n", las=1, cex.points=3, pch=NA, 
#'      col.labels = "red", show.labels=TRUE, cex.0=0.2, 
#'      show.0 = TRUE, col.0="blue", pch.0=4)
#' plot(o_mu1p2_NB, result="dataOCF", category=NA, 
#'      bty="n", las=1)
#' plot(o_mu1p2_NB, result="dataECF", category=NA, 
#'      bty="n", las=1)
#'      
#' plot(o_mu1p2_NB, result="CF", bty="n", las=1)
#' 
#' plot(o_mu1p2_NB, result="OCF", category=1, bty="n", las=1)
#' plot(o_mu1p2_NB, result="OCF", category=2, bty="n", las=1)
#' 
#' plot(o_mu1p2_NB, result="ECFOCF", bty="n", las=1)
#' 
#' plot(o_mu1p2_NB, result="ECFOCF0", bty="n", las=1)
#' plot(o_mu1p2_NB, result="ECFOCF0", category=1, bty="n", las=1)
#' plot(o_mu1p2_NB, result="ECFOCF0", category=2, bty="n", las=1)
#' 
#' plot(o_mu1p2_NB, result="Prob", category=c(1, 2), bty="n", las=1)
#' plot(o_mu1p2_NB, result="Prob", category=c(2, 1), bty="n", las=1)
#' 
#' }
#' @method plot ECFOCF
#' @export


# plot de la table ECF OCF ####

plot.ECFOCF <- function(x, ..., result="CF", category=NA, period=1) {
  p3p <- list(...)
  
  result <- tolower(result)
  
  if (result=="data") {
    do.call(getFromNamespace("plot.TableECFOCF", ns="phenology"), modifyList(list(x=x$data, 
                                  period=period, 
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
                                  show.scale = TRUE), p3p))
  }

  if (result=="ecfocf0") {
    if (all(is.na(category)) | (all(category == ""))) {
      do.call(getFromNamespace("plot.TableECFOCF", ns="phenology"), modifyList(list(x=x$ECFOCF_0, 
                                    period=period, 
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
                                    show.scale = TRUE), p3p))
    } else {
      do.call(getFromNamespace("plot.TableECFOCF", ns="phenology"), 
              modifyList(list(x=x$ECFOCF_0_categories[[as.numeric(category)]], 
                                    period=period, 
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
                              show.scale = TRUE), p3p))
    }
  }
  
  if (result=="ecfocf") {
    if (all(is.na(category)) | (all(category == ""))) {
      do.call(getFromNamespace("plot.TableECFOCF", ns="phenology"), 
              modifyList(list(x=x$ECFOCF, 
                                    period=period, 
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
                              show.scale = TRUE), p3p))
    } else {
      do.call(getFromNamespace("plot.TableECFOCF", ns="phenology"), 
              modifyList(list(x=x$ECFOCF_categories[[as.numeric(category)]],
                              period=period, 
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
                              show.scale = TRUE), p3p))
    }
  }
  
  if (result=="cf") {
    if (all(is.na(category)) | (all(category == ""))) {
        cf <- x$CF
       main="Clutch Frequency: All categories"
    } else {
      cf <- x$CF_categories[[as.numeric(category)]]
      main=paste0("Clutch Frequency: Category ", as.character(category))
    }
    do.call(plot, modifyList(list(x=1:length(cf), 
                                  xlab="Clutch Frequency", 
                                  ylab="Density", 
                                  main=main, 
                                  y=cf, type="h", xaxt="n"), p3p)[c("x", "y", "type", "col", 
                                                          "main", "cex.axis", "bty", "las", 
                                                          "xlab", "ylab", "xaxt", "xlim", "ylim")])
    axis(side = 1, at=1:length(cf), cex.axis=unlist(modifyList(list(cex.axis=0.8), p3p)[c("cex.axis")]))
  }
  
  if ((result=="ocf") | (result=="dataocf")) {
    if (result=="ocf") {
      ylab="Density"
      if (all(is.na(category)) | (all(category == ""))) {
        ecfocf <- x$ECFOCF[, , period]
        main="Observed Clutch Frequency: All categories"
      } else {
        ecfocf <- x$ECFOCF_categories[[as.numeric(category)]][, , period]
        main=paste0("Observed Clutch Frequency: Category ", as.character(category))
      }
    } else {
      ylab="Frequency"
      ecfocf <- x$data[, , period]
      main="Observed OCF"
    }
    ocf <- rowSums(ecfocf, na.rm = TRUE)
    do.call(plot, modifyList(list(x=0:(length(ocf)-1), 
                                  xlab="Observed Clutch Frequency", 
                                  ylab=ylab, 
                                  main=main,
                                  y=ocf, type="h"), p3p)[c("x", "y", "type", "col", 
                                                           "main", "cex.axis", "bty", "las", 
                                                           "xlab", "ylab", "xlim", "ylim", 
                                                           "xaxt", "yaxt", "axes")])
  }
  
  if ((result=="ecf") | (result=="dataecf")) {
    if (result=="ecf") {
      ylab="Density"
      if (all(is.na(category)) | (all(category == ""))) {
        ecfocf <- x$ECFOCF[, , period]
        main="Estimated Clutch Frequency: All categories"
      } else {
        ecfocf <- x$ECFOCF_categories[[as.numeric(category)]][, , period]
        main=paste0("Estimated Clutch Frequency: Category ", as.character(category))
      }
    } else {
      ylab="Frequency"
      ecfocf <- x$data[, , period]
      main="Observed ECF"
    }
    ecf <- colSums(ecfocf, na.rm = TRUE)
    do.call(plot, modifyList(list(x=0:(length(ecf)-1), 
                                  xlab="Estimated Clutch Frequency", 
                                  ylab=ylab, 
                                  main=main, 
                                  y=ecf, type="h"), p3p)[c("x", "y", "type", "col", 
                                                           "main", "cex.axis", "bty", "las", 
                                                           "xlab", "ylab", "xlim", "ylim", 
                                                           "xaxt", "yaxt", "axes")])
    
  }
  
  if (result=="period") {
    if (all(is.na(category)) | (all(category == ""))) {
      y <- x$period
      main="All categories"
    } else {
      y <- x$period_categories[[as.numeric(category)]]
      main=paste("Category", category)
    }
    
    perr <- list(x=0:(length(y)-1), 
                 y=y, 
                 las=1, bty="n", 
                 ylab="Probability of nesting", xlab="Period", 
                 main=main)
    perr <- modifyList(perr, p3p)
    
    do.call(plot, perr)
  }
  
  if (result=="prob") {
    if (is.null(x$SE_df)) {
      warning("The estimate of standard error for capture probability is not available")
    } else {
    if (all(is.na(category)) | (all(category == ""))) {
     category <- ""
     
     cl <- sapply(X = rownames(x$SE_df), function(x) {
       # Je le prends si probx
       grepl("prob", x)
       })
     
    } else {
      # Si j'ai une catégorie
    category <- as.character(category)
    
    cl <- sapply(X = rownames(x$SE_df), function(x) {
      # Je le prends si probx
      grepl(paste0("prob", category, "\\."), x) | 
        # ou ax
        grepl(paste0("a", category), x) | 
        # ou prob. et category = 1 mais pas si ax
        ((grepl(paste0("prob\\."), x)) & (category == "1") & !grepl(paste0("a[^", category, "]"), x))
      })
}
    
    perr <- list(x=1:sum(cl), 
                 y=x$SE_df[cl, "Estimate"], 
                 y.minus=x$SE_df[cl, "2.5 %"], 
                 y.plus = x$SE_df[cl, "97.5 %"], 
                 xlim=c(0.5, sum(cl)+0.5), ylim=c(0,1), 
                 las=1, bty="n", xaxt="n", 
                 ylab="Probability of capture", xlab="Categories", 
                 main="SE using delta method")
    perr <- modifyList(perr, p3p)
    
    do.call(plot_errbar, perr)
    
    if (is.null(p3p$xaxt)) p3p$xaxt <- "r"
    
    if (p3p$xaxt != "n") {
    segments(x0=1:sum(cl), 
             y0=-0.1, y1=-0.15, xpd=TRUE)
    cex <- p3p[["cex.axis"]]
    y <- p3p[["y.axis"]]
    if (is.null(cex)) cex <- 1
    if (is.null(y)) y <- -0.3
    do.call(text, modifyList(list(x = 1:sum(cl), 
         y=y, 
         cex=cex, 
         labels = rownames(x$SE_df)[cl], 
         xpd=TRUE), p3p[c("srt", "labels")]))
    }
    }
  }
}
