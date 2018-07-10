#' IPModel estimates the pattern of internesting intervals for a set of parameters.
#' @title Estimates the pattern of internesting intervals for a set of parameters.
#' @author Marc Girondot
#' @return Return a list with two elements.\cr
#' @param par Set of parameters
#' @param parallel If TRUE, will use parallel computing
#' @param limits A list of limits for various parameters
#' @description This function fits a model of internesting period.\cr
#' The parameters are:\cr
#' \itemize{
#'   \item \code{meanIP} : The average number of days between two nesting processes
#'   \item \code{DeltameanIP} : The shift in days for IP at each new clutch.
#'   \item \code{sdIP} : The standard deviation of number of days between two nesting processes
#'   \item \code{minIP} : The minimum number of days between two nesting processes
#'   \item \code{pAbort} : The -logit of the probability to abort a nesting process
#'   \item \code{meanAbort} : The average of the number of days after the abortion of a nesting process
#'   \item \code{sdAbort} : The standard deviation of the number of days after the abortion of a nesting process
#'   \item \code{pCapture} : The -logit of the probability to capture a female on the beach
#'   \item \code{meanECF} : The average number of clutch a female will try to do being reprensented as ECF
#'   \item \code{sdECF} : The standard deviation of number of clutch a female will try to do
#'   \item \code{N} : The number of replicates to generate the distribution (default is 10000 if not indicated)
#'   \item \code{ECF.x} : The relative proportion of females nesting with ECF = x (ECF.1 being fixed to 1)
#' }
#'  
#' @family Model of Internesting Period
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' par <- c(meanIP = 9.8, 
#' sdIP = 0.1, 
#' minIP = 7, 
#' 
#' pAbort = -logit(0.1), 
#' meanAbort = 2, 
#' sdAbort = 0.05, 
#' 
#' pCapture = -logit(0.8), 
#' 
#' meanECF = 4, 
#' sdECF = 0.1)
#' 
#' model <- IPModel(c(par, N=10000))
#' 
#' plot(model)
#' 
#' }
#' @export

# Lancement du fit ####

IPModel <- function (par, parallel = TRUE, limits = list(meanIP = 40, meanECF = 4, 
                                                         minIP = 15, sdAbort = 1, sdIP = 1, sdECF = 1, DeltameanIP = 0.5, 
                                                         maxDays = 365)) 
{
  meanIP <- abs(par["meanIP"])
  sdIP <- abs(par["sdIP"])
  minIP <- abs(par["minIP"])
  DeltameanIP <- par["DeltameanIP"]
  if (is.na(DeltameanIP)) 
    DeltameanIP <- 0
  pAbort <- invlogit(-par["pAbort"])
  meanAbort <- log(abs(par["meanAbort"]))
  sdAbort <- abs(par["sdAbort"])
  pCapture <- invlogit(-par["pCapture"])
  meanECF <- log(abs(par["meanECF"]))
  sdECF <- abs(par["sdECF"])
  ECF <- abs(par[substr(names(par), 1, 3) == "ECF"])
  N <- par["N"]
  if (is.na(N)) 
    N <- 1e+06
  if (is.na(sdECF)) {
    ECF["ECF.1"] <- 1
    lmx <- ECF[order(as.numeric(gsub("ECF\\.", "", names(ECF))))]
    lmx <- lmx/sum(lmx)
    lmx <- cumsum(lmx)
    meanECF <- log(length(lmx))
    sdECF <- 0
  }
  nc <- 1
  if ((parallel) & (.Platform$OS.type == "unix")) 
    nc <- detectCores()
  ind <- rep(N/nc, nc - 1)
  ind <- c(ind, N - sum(ind))
  if ((meanIP < limits$meanIP) & (meanECF < limits$meanECF) & 
      (minIP < limits$minIP) & (sdAbort < limits$sdAbort) & 
      (sdIP < limits$sdIP) & (sdECF < limits$sdECF) & (abs(DeltameanIP) < 
                                                       limits$DeltameanIP)) {
    rmcl <- mclapply(X = ind, mc.cores = nc, FUN = function(Nx) {
      cumuld <- rep(0, limits$maxDays + 1)
      if (sdECF != 0) {
        CFx <- floor(rlnorm(Nx, meanlog = meanECF, sdlog = sdECF)) +1
      } else {
        CFx <- findInterval(runif(Nx), lmx)+1
      }
      reverseECF <- matrix(data = 0, nrow = 20, ncol = limits$maxDays)
      
      for (individual in 1:Nx) {
        
        CF <- CFx[individual]
        
        d <- 1
        premier <- TRUE
        repeat {
          if (runif(1) <= pAbort) {
            di <- floor(rlnorm(1, meanlog = meanAbort, 
                               sdlog = sdAbort))
            d <- d + di
            if ((runif(1) <= pCapture) & (d < limits$maxDays + 1)) {
              cumuld[d] <- cumuld[d] + 1
              reverseECF[1, d] <- reverseECF[1, d] + 1
            }
            premier <- FALSE
          } else {
            if ((runif(1) <= pCapture) & (d < limits$maxDays + 1) & (!premier)) {
              cumuld[d] <- cumuld[d] + 1
              reverseECF[1, d] <- reverseECF[1, d] +  1
            }
            break
          }
        }
        if (CF > 1) {
          for (CFi in 2:CF) {
            mp <- meanIP + DeltameanIP * (CFi - 2)
            if (mp <= 0)  break
            cpt <- 1
            repeat {
              di <- floor(rlnorm(1, meanlog = log(mp), sdlog = sdIP))
              cpt <- cpt + 1
              if ((di >= minIP) | (cpt > 10)) 
                break
            }
            if (cpt > 10) break
            d <- d + di
            repeat {
              if (runif(1) <= pAbort) {
                di <- floor(rlnorm(1, meanlog = meanAbort, 
                                   sdlog = sdAbort))
                d <- d + di
                if ((runif(1) <= pCapture) & (d < limits$maxDays + 
                                              1)) {
                  cumuld[d] <- cumuld[d] + 1
                  if (CFi < 21) 
                    reverseECF[CFi, d] <- reverseECF[CFi, d] + 1
                }
              } else {
                if ((runif(1) <= pCapture) & (d < limits$maxDays +  1)) {
                  cumuld[d] <- cumuld[d] + 1
                  if (CFi < 21) 
                    reverseECF[CFi, d] <- reverseECF[CFi,  d] + 1
                }
                break
              }
            }
          }
        }
        if (d > limits$maxDays + 1) 
          break
      }
      return(list(cumuld = cumuld, reverseECF = reverseECF))
    })
    cumuld <- rep(0, limits$maxDays + 1)
    reverseECF <- matrix(data = 0, nrow = 20, ncol = limits$maxDays)
    for (j in 1:length(rmcl)) {
      cumuld <- cumuld + rmcl[[j]]$cumuld
      reverseECF <- reverseECF + rmcl[[j]]$reverseECF
    }
  } else {
    cumuld <- rep(1, limits$maxDays + 1)
    reverseECF <- matrix(data = 1, nrow = 20, ncol = limits$maxDays)
  }
  maxd <- which(cumuld != 0)
  if (length(maxd) != 0) {
    maxd <- max(maxd)
    cumuld <- cumuld[1:maxd]
  }
  cumuld <- cumuld/sum(cumuld)
  if (length(cumuld) == 0) {
    cumuld <- 1
    names(cumuld) <- as.character(0:(length(cumuld) - 1))
    reverseECF <- NA
  } else {
  names(cumuld) <- as.character(0:(length(cumuld) - 1))
  colnames(reverseECF) <- paste0("Day", 0:(limits$maxDays - 
                                             1))
  rownames(reverseECF) <- paste0("Cluth", 0:19)
  for (col in 1:limits$maxDays) reverseECF[, col] <- reverseECF[, 
                                                                col]/sum(reverseECF[, col])
  reverseECF[is.na(reverseECF)] <- 0
  
  }
  out <- list(cumuld = cumuld, reverseECF = reverseECF)
  class(out) <- "IP"
  return(out)
}
