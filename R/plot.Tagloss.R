#' plot.tagloss plots the daily rate of tag loss.
#' @title Plot the daily rate of tag loss.
#' @author Marc Girondot
#' @return A dataframe with values used for plotting.
#' @param x Object obteined from Tagloss_fit()
#' @param t Time for which values of model must be ploted
#' @param fitted.parameters Set of parameters
#' @param fixed.parameters Another set of parameters without standard error associated
#' @param model_before Transformation of parameters before to use Tagloss_model()
#' @param model_after Transformation of parameters after to use Tagloss_model()
#' @param model Can be 1, 2, R1, R2, L1, L2 or Cumul (2 tags) or Cumul1 (1 tag)
#' @param col The colors of shading areas of cumul or the color of line
#' @param text.col The text color for cumul model
#' @param scale Scale value. When Cumul is used, scale is always 1.
#' @param add Should the data be added to a previous plot?
#' @param hessian Should confidence interval be shown
#' @param replicates Number of replicates for confidence interval
#' @param probs Quantiles to show for confidence interval
#' @param decoration Try to add name of parameters on the graph
#' @param progressbar Is shown a progressbar?
#' @param ... Parameters transmitted to plot
#' @description Plot the daily rate of tag loss.\cr
#' To use this function without a result of Tagloss_fit(), see the hack in examples.\cr
#' It does not work still with cumul and LR model.
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
#' phenology:::plot.Tagloss(x=NULL, t=t, fitted.parameters=par, model="1", 
#'                          scale=1000, decoration = TRUE)
#' 
#' par <- c(D1_2=200, D2D1_2=100, D3D2_2=200, 
#'          A_2=-logit(0.05), B_2=-logit(0.03), C_2=-logit(0.03))
#' phenology:::plot.Tagloss(x=list(), t=t, fitted.parameters=par, ylim=c(0, 1), 
#'                          scale = 10, model="2")
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
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, model="Cumul", 
#'                          decoration = TRUE)
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, model="Cumul", 
#'                          decoration = TRUE, col=c("red", "green", "blue"))
#'                          
#' # To plot the history of individuals
#' par <- c(D1_R1=200, D2D1_R1=100, D3D2_R1=200, 
#'          A_R1=-logit(5E-4), B_R1=-logit(4E-4), C_R1=-logit(5E-4), 
#'          D1_R2=200, D2D1_R2=100, D3D2_R2=200, 
#'          A_R2=-logit(6E-4), B_R2=-logit(5E-4), C_R2=-logit(6E-4), 
#'          D1_L1=200, D2D1_L1=100, D3D2_L1=200, 
#'          A_L1=-logit(5E-4), B_L1=-logit(4E-4), C_L1=-logit(5E-4), 
#'          D1_L2=200, D2D1_L2=100, D3D2_L2=200, 
#'          A_L2=-logit(6E-4), B_L2=-logit(5E-4), C_L2=-logit(6E-4))
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, model="Cumul", 
#'                          decoration = TRUE)

#' }
#' @method plot Tagloss
#' @export


plot.Tagloss <- function(x, t=NULL, fitted.parameters=NULL, fixed.parameters=NULL, scale = 1, 
                         model_before=NULL,
                         model_after=NULL,
                         model=c("1", "2", "R1", "R2", "L1", "L2", "cumul", "cumul1"), 
                         col=grey.colors(4), 
                         text.col=rev(grey.colors(4)), 
                         add=FALSE, 
                        hessian=NULL, replicates=1000, 
                        probs=c(0.025, 0.975), 
                        progressbar=FALSE, 
                        decoration=FALSE, 
                        ...) {
  
  p3p <- list(...)
  
  # p3p <- list(); x <- NULL; t=NULL; fitted.parameters=NULL; fixed.parameters=NULL; scale = 1; model_before=NULL;model_after=NULL;model="1"; add=FALSE; hessian=NULL; replicates=1000; probs=c(0.025, 0.975); progressbar=FALSE; decoration=FALSE
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
      
      if (!is.null(model_before)) eval(parse(text=model_before), envir= environment())
      
      if (any(grepl("_R2", names(par)))) pR2 <- Tagloss_model(t, par[grepl("_R2", names(par))]) else pR2 <- NA
      if (any(grepl("_R1", names(par)))) pR1 <- Tagloss_model(t, par[grepl("_R1", names(par))]) else pR1 <- NA
      if (any(grepl("_L2", names(par)))) pL2 <- Tagloss_model(t, par[grepl("_L2", names(par))]) else pL2 <- NA
      if (any(grepl("_L1", names(par)))) pL1 <- Tagloss_model(t, par[grepl("_L1", names(par))]) else pL1 <- NA
      
      if (any(grepl("_2", names(par)))) p2 <- Tagloss_model(t, par[grepl("_2", names(par))]) else p2 <- NA
      if (any(grepl("_1", names(par)))) p1 <- Tagloss_model(t, par[grepl("_1", names(par))]) else p1 <- NA
      
      if (!is.null(model_after)) eval(parse(text=model_after), envir= environment())
      
      if (is.na(p1[1]) & !is.na(p2[1])) p1 <- p2
      if (is.na(p2[1]) & !is.na(p1[1])) p2 <- p1
      
      if (is.na(pR1[1]) & !is.na(pR2[1])) pR1 <- pR2
      if (is.na(pR2[1]) & !is.na(pR1[1])) pR2 <- pR1
      
      if (is.na(pL1[1]) & !is.na(pL2[1])) pL1 <- pL2
      if (is.na(pL2[1]) & !is.na(pL1[1])) pL2 <- pL1
      
      
      if ((!is.na(p1[1])) & (!is.na(p2[1]))) {
        d <- data.frame(time=t, N2=rep(x = 1, times=days.maximum), N1=rep(x=0, times=days.maximum), N0=rep(x=0, times=days.maximum))
        for (i in 2:days.maximum) {
          d[i, "N2"] <- d[i-1, "N2"] * (1-p2[i-1])
          d[i, "N1"] <- d[i-1, "N1"] * (1-p1[i-1]) + d[i-1, "N2"] * p2[i-1] * (1-p1[i-1])
          d[i, "N0"] <- d[i-1, "N0"] + d[i-1, "N1"] * p1[i-1] + d[i-1, "N2"] * p2[i-1] * p1[i-1]
        }
        d[, -1] <- d[, -1]/rowSums(d[, -1])
        
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
        d <- data.frame(time=t, NLR=rep(x = 1, times=days.maximum), NL0=rep(x=0, times=days.maximum), N0R=rep(x=0, times=days.maximum), N00=rep(x=0, times=days.maximum))
        for (i in 2:days.maximum) {
          d[i, "NLR"] <- d[i-1, "NLR"] * (1-pR2[i-1]) * (1-pL2[i-1])
          d[i, "NL0"] <- d[i-1, "NL0"] * (1-pR1[i-1]) + d[i-1, "NLR"] * pL2[i-1] * (1-pR1[i-1])
          d[i, "N0R"] <- d[i-1, "N0R"] * (1-pL1[i-1]) + d[i-1, "NLR"] * pR2[i-1] * (1-pL1[i-1])
          d[i, "N00"] <- d[i-1, "N00"] + d[i-1, "NL0"] * pL1[i-1] + d[i-1, "N0R"] * pR1[i-1] + d[i-1, "NLR"] * pL2[i-1] * pR1[i-1] + d[i-1, "NLR"] * pR2[i-1] * pL1[i-1]
        }
        d[, -1] <- d[, -1]/rowSums(d[, -1])
        
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
  # Je plot le modèle des pertes de bagues
  if (progressbar & !is.null(hessian) & replicates>1) {
    if (is.element('progress', installed.packages()[,1])) {
      # library("progress")
      pb <- getFromNamespace("progress_bar", ns="progress")$new(
        format = "  completion [:bar] :percent eta: :eta",
        total = replicates, clear = FALSE)
      libp <- TRUE
    } else {
      libp <- FALSE
      pb <- txtProgressBar(min=0, max=replicates, style=3)
    }
  }

  if (!is.null(hessian)) {
  
  SE <- SEfromHessian(hessian)[names(fitted.parameters)]
  
  setagloss <- matrix(data = rep(x = NA, replicates*length(t)), ncol=replicates)
  # parametre <- matrix(data = rep(x = NA, replicates*length(SE)), ncol=replicates)
  
  # plot(x=t, y=tagloss_f(t, c(par, pfixed)), type="l", bty="n", ylim=c(0, 0.003))
  
  for (replicate in 1:replicates) {
    
    if (progressbar & !is.null(hessian) & replicates>1) {
      if (libp) pb$tick() else setTxtProgressBar(pb, replicate)
    }
    
    
    par <- rnorm(length(SE), fitted.parameters, SE)
    names(par) <- names(SE)
    
    par <- c(par, fixed.parameters)
    
    if (!is.null(model_before)) eval(parse(text=model_before), envir= environment())
    
    if (any(grepl("_R2", names(par)))) pR2 <- Tagloss_model(t, par[grepl("_R2", names(par))]) else pR2 <- NA
    if (any(grepl("_R1", names(par)))) pR1 <- Tagloss_model(t, par[grepl("_R1", names(par))]) else pR1 <- NA
    if (any(grepl("_L2", names(par)))) pL2 <- Tagloss_model(t, par[grepl("_L2", names(par))]) else pL2 <- NA
    if (any(grepl("_L1", names(par)))) pL1 <- Tagloss_model(t, par[grepl("_L1", names(par))]) else pL1 <- NA
    
    if (any(grepl("_2", names(par)))) p2 <- Tagloss_model(t, par[grepl("_2", names(par))]) else p2 <- NA
    if (any(grepl("_1", names(par)))) p1 <- Tagloss_model(t, par[grepl("_1", names(par))]) else p1 <- NA
    
    if (!is.null(model_after)) eval(parse(text=model_after), envir= environment())
    
    if (is.na(p1[1]) & !is.na(p2[1])) p1 <- p2
    if (is.na(p2[1]) & !is.na(p1[1])) p2 <- p1
    
    if (is.na(pR1[1]) & !is.na(pR2[1])) pR1 <- pR2
    if (is.na(pR2[1]) & !is.na(pR1[1])) pR2 <- pR1
    
    if (is.na(pL1[1]) & !is.na(pL2[1])) pL1 <- pL2
    if (is.na(pL2[1]) & !is.na(pL1[1])) pL2 <- pL1
    
    
    if (model == "1") Tlm <- p1
    if (model == "2") Tlm <- p2
    if (model == "r2") Tlm <- pR2
    if (model == "r1") Tlm <- pR1
    if (model == "l2") Tlm <- pL2
    if (model == "l1") Tlm <- pL1
  
    setagloss[, replicate] <- Tlm
  }
  
  mean_ <- apply(setagloss, MARGIN = 1, FUN=mean)
  SE_ <- apply(setagloss, MARGIN = 1, FUN=quantile, probs=probs)
  }
  
  
  par <- c(fitted.parameters, fixed.parameters)
  
  if (!is.null(model_before)) eval(parse(text=model_before), envir= environment())
  
  if (any(grepl("_R2", names(par)))) pR2 <- Tagloss_model(t, par[grepl("_R2", names(par))]) else pR2 <- NA
  if (any(grepl("_R1", names(par)))) pR1 <- Tagloss_model(t, par[grepl("_R1", names(par))]) else pR1 <- NA
  if (any(grepl("_L2", names(par)))) pL2 <- Tagloss_model(t, par[grepl("_L2", names(par))]) else pL2 <- NA
  if (any(grepl("_L1", names(par)))) pL1 <- Tagloss_model(t, par[grepl("_L1", names(par))]) else pL1 <- NA
  
  if (any(grepl("_2", names(par)))) p2 <- Tagloss_model(t, par[grepl("_2", names(par))]) else p2 <- NA
  if (any(grepl("_1", names(par)))) p1 <- Tagloss_model(t, par[grepl("_1", names(par))]) else p1 <- NA
  
  if (!is.null(model_after)) eval(parse(text=model_after), envir= environment())
  
  if (is.na(p1[1]) & !is.na(p2[1])) p1 <- p2
  if (is.na(p2[1]) & !is.na(p1[1])) p2 <- p1
  
  if (is.na(pR1[1]) & !is.na(pR2[1])) pR1 <- pR2
  if (is.na(pR2[1]) & !is.na(pR1[1])) pR2 <- pR1
  
  if (is.na(pL1[1]) & !is.na(pL2[1])) pL1 <- pL2
  if (is.na(pL2[1]) & !is.na(pL1[1])) pL2 <- pL1
  
  
  if (model == "1") {Tlm <- p1; par <- par[grepl("_1", names(par))]}
  if (model == "2") {Tlm <- p2; par <- par[grepl("_2", names(par))]}
  if (model == "r2") {Tlm <- pR2; par <- par[grepl("_R2", names(par))]}
  if (model == "r1") {Tlm <- pR1; par <- par[grepl("_R1", names(par))]}
  if (model == "l2") {Tlm <- pL2; par <- par[grepl("_L2", names(par))]}
  if (model == "l1") {Tlm <- pL1; par <- par[grepl("_L1", names(par))]}
  
  max_y <- max(Tlm)*scale
  if (!is.null(hessian)) {
    max_y <- max(c(max_y, SE_[1, ], SE_[2, ]))
  }
 
  d <- data.frame(time=t, p=Tlm)
  
  if (add) {
    p3px <- modifyList(list(x=t, y=Tlm*scale, col=col[1]), p3p)
    do.call(lines, p3px)
  } else {
    p3px <- modifyList(list(x=t, y=Tlm*scale, las=1, 
                           type="l", ylim=c(0, max_y), 
                           bty="n", col=col[1], 
                           xlab="Days after tagging", 
                           ylab=paste0("Daily tag loss (x", scale, ")")), p3p)
    do.call(plot, p3px)
  }
  
  if (!is.null(hessian)) {
    d <- cbind(d, lower=SE_[1, ], upper=SE_[2, ])
    p3px <- modifyList(list(x=t, y=(SE_[1, ])*scale, lty=2), p3p)
    do.call(lines, p3px)
    p3px <- modifyList(list(x=t, y=(SE_[2, ])*scale, lty=2), p3p)
    do.call(lines, p3px)
  }
  
  
  if (decoration) {
    Arrows <- getFromNamespace(".Arrows", ns="HelpersMG")
    names(par) <- gsub("_.*", "", names(par))
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
    
    if (is.na(par["delta"])) delta <- 1 else delta <- invlogit(-par["delta"])
    if (D1==0) D1_2 <- 0.1 else D1_2 <- 2*D1
    if (D3==D2) D3D2 <- (D3-0.1)*D2 else D3D2 <-D3-D2
    

    A <- A * delta
    B <- B * delta
    C <- C * delta
    
    if (!is.na(par["A"])) {
      Arrows(x1=D1/10, x0=D1/2, y0=A*scale, y1=A*scale)
      text(D1/2+D1/10, y=A*scale, labels = "A", pos=4, col=text.col[1])
    }
    
    if (!is.na(par["B"])) {
    Arrows(x1=D1-D1/10, x0=D1-D1/2, y0=B*scale, y1=B*scale)
    if (B!=C) {
      text(D1-D1/2-D1/20, y=B*scale, labels = "B", pos=2, col=text.col[1])
      if (!is.na(par["C"])) {
        Arrows(x1=D3-D3/10, x0=D3-D3/2, y0=C*scale, y1=C*scale)
        text(D3-D3/2-D3/10, y=C*scale, labels = "C", pos=2, col=text.col[1])
      }
    } else {
      text(D1-D1/2-D1/20, y=B*scale, labels = "B & C", pos=2, col=text.col[1])
    }
    }
    
    mx <- ScalePreviousPlot()$ylim["end"]
    
    par(xpd=TRUE)
    lines(x = c(D1, D1), y=c(0, mx), lty=2)
    if (D1 != D2) {
      text(x=D1+max(t)/70, y=0*scale, labels="D1", pos=4, col=text.col[1])
      lines(x = c(D2, D2), y=c(0, mx), lty=2)
      text(x=D2+max(t)/70, y=0*scale, labels="D2", pos=4, col=text.col[1])
    } else {
      text(x=D1+max(t)/70, y=0*scale, labels="D1 & D2", pos=4, col=text.col[1])
    }
    lines(x = c(D3, D3), y=c(0, mx), lty=2)
    text(x=D3+max(t)/70, y=0*scale, labels="D3", pos=4, col=text.col[1])
  }
  }
  return(invisible(d))
}