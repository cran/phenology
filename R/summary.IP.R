#' summary.IP prints the information from a IP object.
#' @title Print the result information from a IP object.
#' @author Marc Girondot
#' @return Nothing
#' @param object A file of class IP
#' @param ... Not used
#' @param N Number of replicates
#' @param probs Probability of confidence interval
#' @description The function summary.IP shows result and estimates confidence interval.
#' @examples
#' \dontrun{
#' library(phenology)
#' # Read a file with data
#' }
#' @method summary IP
#' @family Model of Internesting Period
#' @export

summary.IP <- function(object, ..., N = NULL, probs=c(0.025, 0.975)) {
  
  ECF <- NULL
  minIP <- NULL
  meanAbort <- NULL
  meanIP <- NULL
  
  level <- qnorm(probs[2])
  
  if (is.list(object)) {
    if (!(identical(object$ML, list())) | !(identical(object$MH, list()))) {
      if (!(identical(object$MH, list()))) {
        data <- object$MH$parametersMCMC$control$data
        pari <- c(as.parameters(object$MH), object$MH$parametersMCMC$control$fixed.parameters)
        SE <- summary(object$MH)$statistics[,"Time-series SE"]
      } else {
        data <- object$ML$data
        pari <- c(object$ML$par, object$ML$fixed.parameters)
        SE <- object$ML$SE
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
        if (is.na(N)) N <- 1000000
      }
      
      if (!is.na(meanECF)) {
        CFx <- floor(rlnorm(N, meanlog=meanECF, sdlog=sdECF))+1
        CFx <- as.data.frame(table(CFx), stringsAsFactors=FALSE)
        CFx[, "CFx"] <- as.numeric(CFx[, "CFx"])
        ECF <- data.frame(ECF=1:max(CFx[, "CFx"]), Freq=0)
        ECF[CFx[, "CFx"], "Freq"] <- CFx[, "Freq"]/sum(CFx[, "Freq"])
      } else {
        ECF["ECF.1"] <- 1
        ECF <- ECF[order(as.numeric(gsub("ECF\\.", "", names(ECF))))]
        ECF <- ECF / sum(ECF)
        ECF <- data.frame(ECF=1:length(ECF), Freq=ECF)
      }
      
        if (!is.null(object$model)  & (Nnull)) {
          model <- object$model
        } else {
          model <- IPModel(pari)
        }
        
        reverseECF <- model$reverseECF
        model <- model$cumuld

      di <- floor(rlnorm(N, meanlog=meanIP, sdlog=sdIP))
      di <- di[di >= minIP]
      di <- as.data.frame(table(di), stringsAsFactors=FALSE)
      di[, "di"] <- as.numeric(di[, "di"])
      IP <- data.frame(IP=1:max(di[, "di"]), Freq=0)
      IP[di[, "di"], "Freq"] <- di[, "Freq"]/sum(di[, "Freq"])
      
      di <- floor(rlnorm(N, meanlog=meanAbort, sdlog=sdAbort))
      di <- as.data.frame(table(di), stringsAsFactors=FALSE)
      di[, "di"] <- as.numeric(di[, "di"])
      Abort <- data.frame(Abort=0:(max(di[, "di"])), Freq=0)
      Abort[di[, "di"]+1, "Freq"] <- di[, "Freq"]/sum(di[, "Freq"])
      
    } else {
      data <- object$cumuld
      reverseECF <- object$reverseECF
    }
  } else {
    pari <- NULL
    data <- object
    reverseECF <- NULL
    SE <- NULL
  }
  
  if (!is.null(pari)) {
    print(paste("Probability of capture", invlogit(-pari["pCapture"])))
    if (!is.null(SE)) {
      # Confidence interval at 95%
      print(paste(invlogit(-pari["pCapture"]-level*SE["pCapture"]), "-", 
                  invlogit(-pari["pCapture"]+level*SE["pCapture"])))
    }
  }
  
  if (!is.null(pari)) {
    print(paste("Probability of aborting nesting process while on the beach", invlogit(-pari["pAbort"])))
    if (!is.null(SE)) {
      # Confidence interval at 95%
      print(paste(invlogit(-pari["pAbort"]-level*SE["pAbort"]), "-", 
                  invlogit(-pari["pAbort"]+level*SE["pAbort"])))
    }
  }
  
  if (!is.null(ECF)) {
    print("Probabilities of ECF")
    print(ECF)
  }
  
  if (!is.null(meanAbort)) {
    print("Mean number of days before a new attempt when a nesting attempt is aborted")
    print(paste("Mean number of days", abs(pari["meanAbort"])))
    if (!is.null(SE)) {
      # Confidence interval at 95%
      print(paste(abs(pari["meanAbort"])-level*SE["meanAbort"], "-", 
                  abs(pari["meanAbort"])+level*SE["meanAbort"]))
    }
  }
  
  di <- rlnorm(N, meanlog=meanAbort, sdlog=sdAbort)
  q <- quantile(di, probs = probs, names=FALSE, type=8)
  print(paste("Confidence interval", probs[1], "-", probs[2], ": ", q[1], "-", q[2]))
  
  if (!is.null(meanIP)) {
    if (!is.na(DeltameanIP)) {
      print("Mean number of days between two nesting attempts (clutch 0)")
    } else {
      print("Mean number of days between two nesting attempts")
    }
    print(paste("Mean number of days", abs(pari["meanIP"])))
    if (!is.null(SE)) {
      # Confidence interval at 95%
      print(paste(abs(pari["meanIP"])-level*SE["meanIP"], "-", 
                  abs(pari["meanIP"])+level*SE["meanIP"]))
    }
  }
  if (!is.na(DeltameanIP)) {
  for (Clutch in 0:(nrow(ECF)-1)) {
    print(paste("Clutch", Clutch))
    mp <- exp(meanIP) + DeltameanIP * Clutch
    di <- rlnorm(N, meanlog=log(mp), sdlog=sdIP)
    di <- di[floor(di)>=minIP]
    q <- unname(quantile(di, probs = probs))
    print(paste("Confidence interval", probs[1], "-", probs[2], ": ", q[1], "-", q[2]))
  }
  } else {
    di <- rlnorm(N, meanlog=meanIP, sdlog=sdIP)
    di <- di[floor(di)>=minIP]
    q <- quantile(di, probs = probs, names=FALSE, type=8)
    print(paste("Confidence interval", probs[1], "-", probs[2], ": ", q[1], "-", q[2]))
  }
  
  

  if (!is.null(minIP)) {
    print("Minimal number of days between two nesting attempts")
    print(paste("Minimal number of days", minIP))
    if (!is.null(SE)) {
      # Confidence interval at 95%
      print(paste(abs(minIP)-level*SE["minIP"], "-", 
                      abs(minIP)+level*SE["minIP"]))
    }
  }
  
  
}
