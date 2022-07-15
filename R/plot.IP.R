#' plot.IP plots a result of Internesting Period fit or data
#' @title Plot a result of Internesting Period fit or data.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Nothing
#' @param x A result for IPFit() or IPModel().
#' @param ... Graphic parameters, see par().
#' @param N Number of replicates for IPModel().
#' @param clutch The rank of clutch when DeltameanIP is used.
#' @param result What result will be plotted: data, model, data&model, IP, Abort, ECF, reverseECF
#' @description This function plots the result of IPFit() or IPModel().\cr
#' If col is defined with a number of colors, only these colors and shown in legend.
#' @family Model of Internesting Period
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' data <- c(0, 47, 15, 6, 5, 4, 2, 5, 57, 203, 205, 103, 35, 24, 12, 10, 
#'   13, 49, 86, 107, 111, 73, 47, 30, 19, 17, 33, 48, 77, 83, 65, 
#'   37, 27, 23, 24, 22, 41, 42, 44, 33, 39, 24, 18, 18, 22, 22, 19, 
#'   24, 28, 17, 18, 19, 17, 4, 12, 9, 6, 11, 7, 11, 12, 5, 4, 6, 
#'   11, 5, 6, 7, 3, 2, 1, 3, 2, 1, 2, 0, 0, 3, 1, 0, 2, 0, 0, 1)
#'   class(data) <- unique(append("IP", class(data)))
#'   plot(data)
#'   
#' ######### Fit parametric ECF using Maximum-Likelihood
#' 
#' par <- c(meanIP = 9.9959691992722917, 
#'          sdIP = 0.10066664270893474, 
#'          minIP = 7.5684588178888754, 
#'          pAbort = 2.2510012544630911, 
#'          meanAbort = 2.8969185085603386, 
#'          sdAbort = 0.92688983853803242, 
#'          pCapture = -1.0393803705929086, 
#'          meanECF = 3.9551519427394255, 
#'          sdECF = 0.31657679943365019)
#' 
#' fML <- IPFit(x=par, 
#' fixed.parameters=c(N=1000000),
#' data=data, 
#' verbose=FALSE, 
#' model="ML")
#' 
#' # Plot the fitted ECF
#' plot(fML, result="ECF")
#' 
#' # Plot the Internesting Period distribution
#' plot(fML, result="IP")
#' 
#' # Plot the distribution of days between tentatives
#' plot(fML, result="Abort")
#' plot(fML, result="Abort", xlim=c(0, 10))
#' 
#' # Plot the data
#' plot(fML, result="data")
#' 
#' # Plot the data and the model
#' plot(fML, result="data&model")
#' 
#' # Plot the cumulative proportion of ECF according to date of observation
#' plot(fML, result="reverseECF")
#' plot(fML_NP_Delta, result="reverseECF", col=grey.colors(128))
#' 
#' ######### Fit using Metropolis-Hastings
#' # ECF.1 = 1 is fixed
#' par <- c(ECF.2 = 0.044151921569961131, 
#'          ECF.3 = 2.0020778325280748, 
#'          ECF.4 = 2.6128345101617083, 
#'          ECF.5 = 2.6450582416622375, 
#'          ECF.6 = 2.715198206774927, 
#'          ECF.7 = 2.0288031327239904, 
#'          ECF.8 = 1.0028041546528881, 
#'          ECF.9 = 0.70977432157689235, 
#'          ECF.10 = 0.086052204035003091, 
#'          ECF.11 = 0.011400419961702518, 
#'          ECF.12 = 0.001825219438794076, 
#'          ECF.13 = 0.00029398731859899116, 
#'          ECF.14 = 0.002784886479846703, 
#'          meanIP = 9.9887100433529721, 
#'          sdIP = 0.10580250625108811, 
#'          minIP = 6.5159124624132048, 
#'          pAbort = 2.5702251748938956, 
#'          meanAbort = 2.2721679285648841, 
#'          sdAbort = 0.52006431730489933, 
#'          pCapture = 0.079471782729506113)
#'          
#' df <- data.frame(Density=rep("dunif", length(par)), 
#' Prior1=c(rep(0, 13), 8, 0.001, 0, -8, 0, 0.001, -8), 
#' Prior2=c(rep(10, 13), 12, 1, 10, 8, 2, 1, 8), 
#' SDProp=unname(c(ECF.2 = 6.366805760909012e-05, 
#'                 ECF.3 = 6.366805760909012e-05, 
#'                 ECF.4 = 6.366805760909012e-05, 
#'                 ECF.5 = 6.366805760909012e-05, 
#'                 ECF.6 = 6.366805760909012e-05, 
#'                 ECF.7 = 6.366805760909012e-05, 
#'                 ECF.8 = 6.366805760909012e-05, 
#'                 ECF.9 = 6.366805760909012e-05, 
#'                 ECF.10 = 6.366805760909012e-05, 
#'                 ECF.11 = 6.366805760909012e-05, 
#'                 ECF.12 = 6.366805760909012e-05, 
#'                 ECF.13 = 6.366805760909012e-05, 
#'                 ECF.14 = 6.366805760909012e-05, 
#'                 meanIP = 6.366805760909012e-05, 
#'                 sdIP = 6.366805760909012e-05, 
#'                 minIP = 6.366805760909012e-05, 
#'                 pAbort = 6.366805760909012e-05, 
#'                 meanAbort = 6.366805760909012e-05, 
#'                 sdAbort = 6.366805760909012e-05, 
#'                 pCapture = 6.366805760909012e-05)),               
#' Min=c(rep(0, 13), 8, 0.001, 0, -8, 0, 0.001, -8), 
#' Max=c(rep(10, 13), 12, 1, 10, 8, 2, 1, 8),
#' Init=par, stringsAsFactors = FALSE)
#' rownames(df)<- names(par)
#' 
#' fMH <- IPFit(parametersMH=df, 
#' fixed.parameters=c(N=10000),
#' data=data, 
#' verbose=FALSE, 
#' n.iter = 10000,
#' n.chains = 1, n.adapt = 100, thin = 1, trace = TRUE,
#' adaptive = TRUE, 
#' model="MH")
#' 
#' # Plot the fitted ECF
#' plot(fMH, result="ECF")
#' 
#' # Plot the posteriors and priors
#' plot(fMH$MH, parameters="meanIP", xlim=c(6, 14))
#' 
#' plot(x=1:length(fMH$MH$resultLnL[[1]]), y=fMH$MH$resultLnL[[1]], 
#' type="l", xlab="Iterations", ylab="Ln L", bty="n", las=1)
#' }
#' @method plot IP
#' @export


plot.IP <- plot.IP <- function (x, ..., N = NULL, clutch = 1, result = "data") {
  
  
  result <- match.arg(arg=result, choices = c("data", "model", "data&model", "IP", "Abort", "ECF", "reverseECF"))
  
  p3p <- list(...)
  result <- tolower(result)
  if (is.list(x)) {
    if (!(identical(x$ML, list())) | !(identical(x$MH, list()))) {
      if (!(identical(x$MH, list()))) {
        data <- x$MH$parametersMCMC$control$data
        pari <- c(as.parameters(x$MH), x$MH$parametersMCMC$control$fixed.parameters)
      } else {
        data <- x$ML$data
        pari <- c(x$ML$par, x$ML$fixed.parameters)
      }
      meanIP <- log(abs(pari["meanIP"]))
      sdIP <- abs(pari["sdIP"])
      minIP <- abs(pari["minIP"])
      DeltameanIP <- pari["DeltameanIP"]
      pAbort <- invlogit(-pari["pAbort"])
      meanAbort <- log(abs(pari["meanAbort"]))
      sdAbort <- abs(pari["sdAbort"])
      pCapture <- invlogit(-pari["pCapture"])
      meanECF <- log(abs(pari["meanECF"]))
      sdECF <- abs(pari["sdECF"])
      ECF <- abs(pari[substr(names(pari), 1, 3) == "ECF"])
      Nnull <- is.null(N)
      if (Nnull) {
        N <- pari["N"]
        if (is.na(N)) 
          N <- 1e+06
      }
      if (!is.na(meanECF)) {
        CFx <- floor(rlnorm(N, meanlog = meanECF, sdlog = sdECF)) + 
          1
        CFx <- as.data.frame(table(CFx), stringsAsFactors = FALSE)
        CFx[, "CFx"] <- as.numeric(CFx[, "CFx"])
        ECF <- data.frame(ECF = 1:max(CFx[, "CFx"]), 
                          Freq = 0)
        ECF[CFx[, "CFx"], "Freq"] <- CFx[, "Freq"]/sum(CFx[, 
                                                           "Freq"])
      } else {
        ECF["ECF.1"] <- 1
        ECF <- ECF[order(as.numeric(gsub("ECF\\.", "", 
                                         names(ECF))))]
        ECF <- ECF/sum(ECF)
        ECF <- data.frame(ECF = 1:length(ECF), Freq = ECF)
      }
      if ((result == "model") | (result == "data&model") | (result == "reverseecf")) {
        if (!is.null(x$model) & (Nnull)) {
          model <- x$model
        } else {
          model <- IPModel(pari)
        }
        reverseECF <- model$reverseECF
        model <- model$cumuld
      } else {
        model <- NULL
        reverseECF <- NULL
      }
      di <- floor(rlnorm(N, meanlog = log(exp(meanIP) + 
                                            (clutch - 1) * ifelse(is.na(DeltameanIP), 0, 
                                                                  DeltameanIP)), sdlog = sdIP))
      
      di <- di[di >= minIP]
      di <- as.data.frame(table(di), stringsAsFactors = FALSE)
      di[, "di"] <- as.numeric(di[, "di"])
      IP <- data.frame(IP = 1:max(di[, "di"]), Freq = 0)
      IP[di[, "di"], "Freq"] <- di[, "Freq"]/sum(di[, "Freq"])
      di <- floor(rlnorm(N, meanlog = meanAbort, sdlog = sdAbort))
      di <- as.data.frame(table(di), stringsAsFactors = FALSE)
      di[, "di"] <- as.numeric(di[, "di"])
      Abort <- data.frame(Abort = 0:(max(di[, "di"])), 
                          Freq = 0)
      Abort[di[, "di"] + 1, "Freq"] <- di[, "Freq"]/sum(di[, 
                                                           "Freq"])
    } else {
      data <- x$cumuld
      reverseECF <- x$reverseECF
    }
  } else {
    pari <- NULL
    data <- x
    reverseECF <- NULL
  }
  if (is.matrix(data)) {
    data <- colSums(data)
  }
  if (result == "data") {
    do.call(plot, modifyList(list(x = (1:length(data)) - 
                                    1, y = data, type = "h", bty = "n", las = 1, xlab = "Days after first observation", 
                                  ylab = "Number of observations"), p3p))
  }
  if (result == "model") {
    do.call(plot, modifyList(list(x = as.numeric(names(model)), 
                                  y = model, type = "h", bty = "n", las = 1, xlab = "Days after first observation", 
                                  ylab = "Number of observations"), p3p))
  }
  if (result == "ecf") {
    do.call(plot, modifyList(list(x = ECF$ECF, y = ECF$Freq, 
                                  type = "h", bty = "n", las = 1, xlab = "Estimated Clutch Frequency", 
                                  ylab = "Proportion", xaxt = "n"), p3p))
    axis(1, 1:(ScalePreviousPlot()$xlim["end"]), cex.axis = 0.8)
  }
  if (result == "data&model") {
    xl <- max(c(as.numeric(names(model)), length(data) - 
                  1))
    if (!is.null(p3p$xlim)) 
      xl <- p3p$xlim[2]
    do.call(plot, modifyList(list(x = (1:length(data)) - 
                                    1, y = data, type = "h", bty = "n", las = 1, ylim = c(-max(c(data, 
                                                                                                 model * sum(data))), +max(c(data, model * sum(data)))), 
                                  xlim = c(0, xl), xaxt = "n", xlab = "Days after first observation", 
                                  ylab = "Number of observations"), p3p))
    par(new = TRUE)
    do.call(plot, modifyList(list(x = as.numeric(names(model)), 
                                  y = -model * sum(data), type = "h", bty = "n", ylim = c(-max(c(data, 
                                                                                                 model * sum(data))), +max(c(data, model * sum(data)))), 
                                  xlim = c(0, xl), axes = FALSE, xlab = "", ylab = "", 
                                  xaxt = "n"), p3p))
    axis(1, 0:(ScalePreviousPlot()$xlim["end"]), cex.axis = 0.8)
    text(xl, max(c(data, model * sum(data)))/2, labels = "Observations", 
         pos = 2)
    text(xl, -max(c(data, model * sum(data)))/2, labels = "Model", 
         pos = 2)
  }
  if (result == "ip") {
    do.call(plot, modifyList(list(x = IP$IP, y = IP$Freq, 
                                  type = "h", bty = "n", las = 1, xlab = "Internesting Period", 
                                  ylab = "Proportion", xaxt = "n"), p3p))
    axis(1, 1:(ScalePreviousPlot()$xlim["end"]), cex.axis = 0.8)
  }
  if (result == "abort") {
    do.call(plot, modifyList(list(x = Abort$Abort, y = Abort$Freq, 
                                  type = "h", bty = "n", las = 1, xlab = "Number of days between two tentatives", 
                                  ylab = "Proportion", xaxt = "n"), p3p))
    axis(1, 0:(ScalePreviousPlot()$xlim["end"]), cex.axis = 0.8)
  }
  if ((result == "reverseecf") & (!is.null(reverseECF))) {
    cmrev <- reverseECF
    maxx <- 1
    maxy <- 1
    for (col in 1:ncol(cmrev)) {
      maxx <- ifelse(test = (sum(cmrev[, col]) != 0), col, 
                     maxx)
      maxy <- max(which(cmrev[, col] != 0), maxy)
      cmrev[, col] <- cumsum(cmrev[, col])
    }
    if (is.null(p3p$col)) {
      color <- rainbow(maxy)
      k_maxy <- maxy
    } else {
      color <- p3p$col
      k_maxy <- length(color)
      if (length(color) < maxy) {
        color <- c(p3p$col, rep(tail(p3p$col, n=1), maxy-length(p3p$col)))
        
        # p3p$col <- color[floor(length(color) * (1:maxy)/(maxy))]
        
      }
    }
    par(mar = c(4, 4, 2, 6) + 0.4)
    par(xpd=FALSE)
    do.call(plot, modifyList(list(x = 0:(ncol(reverseECF) - 
                                           1), y = reverseECF[1, ], bty = "n", ylim = c(0, 1), 
                                  xlim = c(0, maxx - 1), las = 1, type = "n", xlab = "Days after first observation", 
                                  ylab = "Cumulative proportion"), p3p))
    polygon(x = c(0:(ncol(reverseECF) - 1), (ncol(reverseECF) - 
                                               1):0), y = c(rep(0, ncol(reverseECF)), rev(cmrev[1, 
                                                                                                ])), col = color[1], border = NA)
    for (i in 2:maxy) {
      polygon(x = c(0:(ncol(reverseECF) - 1), (ncol(reverseECF) - 
                                                 1):0), y = c(cmrev[i - 1, ], rev(cmrev[i, ])), 
              col = color[i], border = NA)
    }
    par(xpd = TRUE)
    maxx <- ScalePreviousPlot()$xlim[2]
    dx <- maxx/10
    for (ecf in 1:k_maxy) {
      y <- 0.1 + (ecf - 1) * 0.9/(k_maxy)
      dy <- 0.1 + (ecf - 1 + 0.8) * 0.9/(k_maxy)
      if (ecf !=1) {
      polygon(x = c(maxx + dx, (maxx + dx)*1.05, (maxx + dx)*1.05, maxx + dx), y = c(y, y, dy, dy), col = color[ecf], 
              border = NA)
      }
      text(x = (maxx + dx)*1.07, y = mean(c(y, dy)), labels = (ecf - 1))
    }
    text((maxx + dx)*1.025, y = 1*1.05, labels = "Clutch")
  }
}
