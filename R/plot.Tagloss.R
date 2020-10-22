#' plot.tagloss plots the daily rate of tag loss.
#' @title Plot the daily rate of tag loss.
#' @author Marc Girondot
#' @return An invisible dataframe with values used for plotting.
#' @param x Object obteined from Tagloss_fit()
#' @param t Time for which values of model must be ploted
#' @param fitted.parameters Set of parameters
#' @param fixed.parameters Another set of parameters without standard error associated
#' @param model_before Transformation of parameters before to use Tagloss_model()
#' @param model_after Transformation of parameters after to use Tagloss_model()
#' @param model Can be 1, 2, R1, R2, L1, L2 or Cumul (2 tags) or Cumul1 (1 tag)
#' @param col The colors of shading areas of cumul or the color of line
#' @param text.col The text color for cumul model
#' @param label.col The text color used for labels when decoration is true
#' @param scale Scale value. When Cumul is used, scale is always 1.
#' @param add Should the data be added to a previous plot?
#' @param hessian Hessian matrix of parameters
#' @param replicates Number of replicates for confidence interval
#' @param probs Quantiles to show for confidence interval
#' @param decoration Try to add name of parameters on the graph
#' @param progressbar Is shown a progressbar?
#' @param ... Parameters transmitted to plot
#' @description Plot the daily rate of tag loss.\cr
#' To use this function without a result of Tagloss_fit(), see the hack in examples.
#' @family Model of Tag-loss
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' t <- 1:1000
#' 
#' par <- c(D1_1=200, D2D1_1=100, D3D2_1=200, 
#'          A_1=-logit(0.02), B_1=-logit(0.05), C_1=-logit(0.07))
#' phenology:::plot.Tagloss(x=list(), t=t, fitted.parameters=par, model="1")
#' phenology:::plot.Tagloss(x=list(), t=t, fitted.parameters=par, model="1", 
#'                          scale=1000, decoration = TRUE)
#' 
#' par <- c(D1_2=200, D2D1_2=100, D3D2_2=200, 
#'          A_2=-logit(0.05), B_2=-logit(0.03), C_2=-logit(0.03))
#' phenology:::plot.Tagloss(x=list(), t=t, fitted.parameters=par, ylim=c(0, 1), 
#'                          scale = 10, model="2", decoration = TRUE)
#' 
#' par <- c(D1_L2=200, D2D1_L2=100, D3D2_L2=200, 
#'          A_L2=-logit(0.05), B_L2=-logit(0.07), C_L2=-logit(0.07))
#' phenology:::plot.Tagloss(x=list(), t=t, fitted.parameters=par, model="L2")
#' 
#' par <- c(D1_R2=200, D2D1_R2=0, D3D2_R2=700, 
#'          A_R2=-logit(0.02), B_R2=-logit(0.05), C_R2=-logit(0.07))
#' phenology:::plot.Tagloss(x=list(), t=t, fitted.parameters=par, model="R2", 
#'                          col="red", add=TRUE)
#' 
#' par <- c(D1_L1=200, D2D1_L1=2000, D3D2_L1=2000, 
#'         A_L1=-logit(0.05), B_L1=-logit(0.02), C_L1=-logit(0.1))
#' phenology:::plot.Tagloss(x=list(), t=t, fitted.parameters=par, model="L1")
#' 
#' # To plot the history of individuals
#' par <- c(D1_1=200, D2D1_1=100, D3D2_1=200, 
#'          A_1=-logit(5E-4), B_1=-logit(4E-4), C_1=-logit(5E-4), 
#'          D1_2=200, D2D1_2=100, D3D2_2=200, 
#'          A_2=-logit(6E-4), B_2=-logit(5E-4), C_2=-logit(6E-4))
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, 
#'                          model="Cumul", 
#'                          decoration = TRUE)
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, 
#'                          model="Cumul", 
#'                          decoration = TRUE, col=c("red", "green", "blue"))
#'                          
#' # To plot the history of individuals
#' par <- c(D1_R1=200, D2D1_R1=300, D3D2_R1=200, 
#'          A_R1=-logit(5E-4), B_R1=-logit(4E-4), C_R1=-logit(5E-4), 
#'          D1_R2=200, D2D1_R2=200, D3D2_R2=200, 
#'          A_R2=-logit(6E-4), B_R2=-logit(5E-4), C_R2=-logit(6E-4), 
#'          D1_L1=200, D2D1_L1=400, D3D2_L1=200, 
#'          A_L1=-logit(5E-4), B_L1=-logit(4E-4), C_L1=-logit(5E-4), 
#'          D1_L2=200, D2D1_L2=100, D3D2_L2=200, 
#'          A_L2=-logit(6E-4), B_L2=-logit(5E-4), C_L2=-logit(6E-4))
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, model="Cumul", 
#'                          decoration = TRUE)
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, model="R1", 
#'                          decoration = TRUE)                         
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, model="R2", 
#'                          decoration = TRUE)                         
#' }
#' @method plot Tagloss
#' @export


plot.Tagloss <- function(x, t=NULL, fitted.parameters=NULL, fixed.parameters=NULL, scale = 1, 
                         model_before=NULL,
                         model_after=NULL,
                         model=c("1", "2", "R1", "R2", "L1", "L2", "cumul", "cumul1", "N2", "N1", "N0", "NLR", "N0R", "NL0", "N00"), 
                         col=rev(grey.colors(4, start = 0.9, end = 0.3)), 
                         text.col=grey.colors(4, start = 0.9, end = 0.3), 
                         label.col="black",
                         add=FALSE, 
                         hessian=NULL, replicates=NULL, 
                         probs=c(0.025, 0.975), 
                         progressbar=FALSE, 
                         decoration=FALSE, 
                         ...) {
  
  p3p <- list(...)
  
  # p3p <- list(); x <- NULL; t=NULL; fitted.parameters=NULL; fixed.parameters=NULL; scale = 1; model_before=NULL;model_after=NULL;model="1"; col=rev(grey.colors(4, start = 0.9, end = 0.3));text.col=grey.colors(4, start = 0.9, end = 0.3);label.col="black";add=FALSE; hessian=NULL; replicates=NULL; probs=c(0.025, 0.975); progressbar=FALSE; decoration=FALSE
  model <- tolower(model[1])
  
  if (is.null(t)) {
    if (!is.null(x)) {
      t <- 1:x$mx
    } else {
      t <- 1:1000
    }
  }
  
  if (length(t)==1) t <- 1:t
  days.maximum <- max(t)
  
  if (is.null(fitted.parameters)) fitted.parameters <- x$par
  if (is.null(fixed.parameters)) fixed.parameters <- x$fixed.par
  if (is.null(hessian)) hessian <- x$hessian
  if (is.null(model_before)) model_before <- x$model_before
  if (is.null(model_after)) model_after <- x$model_after
  
  if ((model=="cumul") | (model=="cumul1")) {
    
    # Je plot le modèle cumulé - Je ne prends pas en compte la variabilité
    par <- c(fitted.parameters, fixed.parameters)
    
    d <- Tagloss_cumul(t=t, par=par, hessian=NULL, model_before=model_before, 
                  model_after = model_after)
    
      if (all(colnames(d) != "NLR")) {
        
      p3px <- modifyList(list(x=c(1, days.maximum), y=c(1, 1), las=1, 
                              type="l", ylim=c(0, 1), 
                              xlim=c(1, days.maximum), 
                              bty="n", 
                              col="black",
                              xlab="Days after tagging", 
                              ylab="Relative frequency", 
                              lty=2), p3p)
      
      p3px <- modifyList(p3px, list(lty=2, border=NULL))
      do.call(plot, p3px)
      
      yl <- c(rep(0, days.maximum), rev(d[, "N2"]))
      xl <- c(t, rev(t))
      
      p3px <- modifyList(list(x=xl, y=yl, col=col[1]), p3p[c("lwd", "lty", "border")])
      p3px <- modifyList(p3px, list(density=NA))
      do.call(polygon, p3px)
      
      if (model=="cumul1") {
        yl <- c(d[, "N2"], rep(1, days.maximum))
        p3px <- modifyList(list(x=xl, y=yl, col=col[2]), p3p[c("lwd", "lty", "border")])
        p3px <- modifyList(p3px, list(density=NA))
        do.call(polygon, p3px)
      } else {
        yl <- c(d[, "N2"], rev(d[, "N2"]+d[, "N1"]))
        p3px <- modifyList(list(x=xl, y=yl, col=col[2]), p3p[c("lwd", "lty", "border")])
        p3px <- modifyList(p3px, list(density=NA))
        do.call(polygon, p3px)
        
        yl <- c(d[, "N2"]+d[, "N1"], rep(1, days.maximum))
        p3px <- modifyList(list(x=xl, y=yl, col=col[3]), p3p[c("lwd", "lty", "border")])
        p3px <- modifyList(p3px, list(density=NA))
        do.call(polygon, p3px)
      }
      
      if (decoration) {
        if (model=="cumul") {
          text(days.maximum/8, d[days.maximum/8, "N2"]/2, labels="2 tags", pos=4, col=text.col[1])
          text(days.maximum/2, d[days.maximum/2, "N2"] + (d[days.maximum/2, "N1"])/2, labels="1 tag", pos=4, col=text.col[2])
          text(days.maximum*7/8, d[days.maximum*7/8, "N2"] + d[days.maximum*7/8, "N1"] + d[days.maximum*7/8, "N0"]/2, labels="0 tag", pos=4, col=text.col[3])
        } else {
          text(days.maximum/2, d[days.maximum/2, "N2"] + (1 - d[days.maximum/2, "N2"])/2, labels="0 tag", pos=4, col=text.col[1])
          text(days.maximum/2, d[days.maximum/2, "N1"] + (d[days.maximum/2, "N2"] - d[days.maximum/2, "N1"])/2, labels="1 tag", pos=4, col=text.col[2])
        }
        
      }
    } else {
      # Modèle R L 

      p3px <- modifyList(list(x=c(1, days.maximum), y=c(1, 1), las=1, 
                              type="l", ylim=c(0, 1), 
                              xlim=c(1, days.maximum), 
                              bty="n", 
                              col="black",
                              xlab="Days after tagging", 
                              ylab="Relative frequency", 
                              lty=2), p3p)
      
      p3px <- modifyList(p3px, list(lty=2, border=NULL))
      do.call(plot, p3px)
      
      yl <- c(rep(0, days.maximum), rev(d[, "NLR"]))
      xl <- c(t, rev(t))
      
      p3px <- modifyList(list(x=xl, y=yl, col=col[1]), p3p[c("lwd", "lty", "border")])
      p3px <- modifyList(p3px, list(density=NA))
      do.call(polygon, p3px)
      
      yl <- c(d[, "NLR"], rev(d[, "NL0"]+d[, "NLR"]))
      xl <- c(t, rev(t))
      
      p3px <- modifyList(list(x=xl, y=yl, col=col[2]), p3p[c("lwd", "lty", "border")])
      p3px <- modifyList(p3px, list(density=NA))
      do.call(polygon, p3px)
      
      yl <- c(d[, "NL0"]+d[, "NLR"], rev(d[, "N0R"]+d[, "NL0"]+d[, "NLR"]))
      xl <- c(t, rev(t))
      
      p3px <- modifyList(list(x=xl, y=yl, col=col[3]), p3p[c("lwd", "lty", "border")])
      p3px <- modifyList(p3px, list(density=NA))
      do.call(polygon, p3px)
      
      yl <- c(d[, "N0R"]+d[, "NL0"]+d[, "NLR"], rep(1, days.maximum))
      xl <- c(t, rev(t))
      
      p3px <- modifyList(list(x=xl, y=yl, col=col[4]), p3p[c("lwd", "lty", "border")])
      p3px <- modifyList(p3px, list(density=NA))
      do.call(polygon, p3px)
      
      if (decoration) {
        text(days.maximum/8, d[days.maximum/8, "NLR"]*2/5, labels="Left and right tags", pos=4, col=text.col[1])
        text(days.maximum*3/8, d[days.maximum*3/8, "NLR"] + (d[days.maximum*3/8, "NL0"])*2/5, labels="Left tag", pos=4, col=text.col[2])
        text(days.maximum*5/8, d[days.maximum*5/8, "NLR"] + d[days.maximum*5/8, "NL0"] + d[days.maximum*5/8, "N0R"]*2/5, labels="Right tag", pos=4, col=text.col[3])
        text(days.maximum*7/8, d[days.maximum*7/8, "NLR"] + d[days.maximum*7/8, "NL0"] + d[days.maximum*7/8, "N0R"] + d[days.maximum*7/8, "N00"]*2/5, labels="No tag", pos=4, col=text.col[4])
      }
      
      
    }
  } else {
    
    # Maintenant j'ai directement la distribution dans Tagloss_model
    par <- c(fitted.parameters, fixed.parameters)
    
    
    Tlm <- Tagloss_model(t, par=par, hessian=hessian, 
                         model_before = model_before, 
                         model=model, replicates=replicates)
    
    if (!is.null(hessian)) {
      Tlm_y <- Tlm[, 2]
      Tlm_l <- Tlm[, 5]
      Tlm_u <- Tlm[, 6]
      
    } else {
      Tlm_y <- Tlm
      Tlm_l <- NULL
      Tlm_u <- NULL
    }
    
    
    
    max_y <- max(c(Tlm_y, Tlm_l, Tlm_u))*scale
    
    d <- data.frame(time=t, p=Tlm_y)
    
    if (add) {
      p3px <- modifyList(list(x=t, y=Tlm_y*scale, col=col[1]), p3p)
      do.call(lines, p3px)
    } else {
      p3px <- modifyList(list(x=t, y=Tlm_y*scale, las=1, 
                              type="l", ylim=c(0, max_y), 
                              bty="n", col=col[1], 
                              xlab="Days after tagging", 
                              ylab=paste0("Daily tag loss (x", scale, ")")), p3p)
      do.call(plot, p3px)
    }
    
    if (!is.null(hessian)) {
      d <- cbind(d, lower=Tlm_l, upper=Tlm_u)
      p3px <- modifyList(list(x=t, y=Tlm_l*scale, lty=2, col=col[1]), p3p)
      do.call(lines, p3px)
      p3px <- modifyList(list(x=t, y=Tlm_u*scale, lty=2, col=col[1]), p3p)
      do.call(lines, p3px)
    }
    
    
    if (decoration) {
      
      if (!is.null(model)) {
        model <- toupper(model)
        par <- par[grepl(paste0(".+_", model), names(par))]
        names(par) <- gsub(paste0("_", model), "", names(par))
      } else {
        names(par) <- gsub("_.*", "", names(par))
      }
      
      #  Je ne le fais que pour le nouveau modèle
      if (any(names(par) == "D1")) {
      
      Arrows <- getFromNamespace(".Arrows", ns="HelpersMG")

      if (!is.na(par["D1"])) D1 <- abs(par["D1"]) else D1 <- 0
      if (!is.na(par["D2D1"])) D2D1 <- abs(par["D2D1"]) else D2D1 <- max(t)+1
      D2 <- D1 + D2D1
      if (!is.na(par["D3D2"])) D3D2 <- abs(par["D3D2"]) else D3D2 <- max(t)+1
      D3 <- D2 + D3D2
      if (!is.na(par["A"])) A <- invlogit(-par["A"]) else A <- 0
      if (A == 0) A <- 1E-9
      if (A == 1) A <- 1-1E-9
      if (!is.na(par["B"])) B <- invlogit(-par["B"]) else B <- 0
      if (B == 0) B <- 1E-9
      if (B == 1) B <- 1-1E-9
      if (!is.na(par["C"])) C <- invlogit(-par["C"]) else C <- 0
      if (C == 0) C <- 1E-9
      if (C == 1) C <- 1-1E-9
      
      if (is.na(par["delta"])) delta <- 0 else delta <- par["delta"]
      if (D1==0) D1_2 <- 0.1 else D1_2 <- 2*D1
      if (D3==D2) D3D2 <- (D3-0.1)*D2 else D3D2 <-D3-D2
      
      
      A <- A + delta
      B <- B + delta
      C <- C + delta
      
      if (!is.na(par["A"])) {
        Arrows(x1=D1/10, x0=D1/2, y0=A*scale, y1=A*scale)
        text(D1/2+D1/10, y=A*scale, labels = "A", pos=4, col=label.col)
      }
      
      if (!is.na(par["B"])) {
        Arrows(x1=D1-D1/10, x0=D1-D1/2, y0=B*scale, y1=B*scale)
        if (B!=C) {
          text(D1-D1/2-D1/20, y=B*scale, labels = "B", pos=2, col=label.col)
          if (!is.na(par["C"])) {
            Arrows(x1=D3-D3/10, x0=D3-D3/2, y0=C*scale, y1=C*scale)
            text(D3-D3/2-D3/10, y=C*scale, labels = "C", pos=2, col=label.col)
          }
        } else {
          text(D1-D1/2-D1/20, y=B*scale, labels = "B & C", pos=2, col=label.col)
        }
      }
      
      mx <- ScalePreviousPlot()$ylim["end"]
      
      par(xpd=TRUE)
      lines(x = c(D1, D1), y=c(0, mx), lty=2)
      if (D1 != D2) {
        text(x=D1+max(t)/70, y=0*scale, labels="D1", pos=4, col=label.col)
        lines(x = c(D2, D2), y=c(0, mx), lty=2)
        text(x=D2+max(t)/70, y=0*scale, labels="D2", pos=4, col=label.col)
      } else {
        text(x=D1+max(t)/70, y=0*scale, labels="D1 & D2", pos=4, col=label.col)
      }
      lines(x = c(D3, D3), y=c(0, mx), lty=2)
      text(x=D3+max(t)/70, y=0*scale, labels="D3", pos=4, col=label.col)
      }
    }
  }
  return(invisible(d))
}
