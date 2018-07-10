#' IPPredict calculates the possible clutch number based on observed Internesting Period.
#' @title Predict the possible clutch number based on observed Internesting Period.
#' @author Marc Girondot
#' @return A data.frame
#' @param x A result of IPFit().
#' @param par A set of parameters.
#' @param IP A vector of Internesting Period
#' @param N Number of replicates
#' @description This function predicts the possible clutch number 
#' based on observed Internesting Period.\cr
#' @family Model of Internesting Period
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' ######### Fit using Maximum-Likelihood
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
#' IPPredict(par=par, IP=c(10, 80))
#' 
#' }
#' @export


IPPredict <- function (x = NULL, par = NULL, N = NULL, IP = 0:100) 
{
  if (!is.null(x)) {
    if (is.null(x$ML)) {
      par <- c(as.parameters(x$MH), x$MH$parametersMCMC$control$fixed.parameters)
    }        else {
      par <- c(x$ML$par, x$ML$fixed.parameters)
    }
    Nnull <- is.null(N)
    if (Nnull) {
      N <- par["N"]
      if (is.na(N)) 
        N <- 1e+06
    }
    par["N"] <- N
    if (!is.null(x$model) & (Nnull)) {
      model <- x$model
    }        else {
      model <- IPModel(par)
    }
  }    else {
    model <- IPModel(par)
  }
  reverseECF <- model$reverseECF
  maxx <- 1
  maxy <- 1
  for (col in 1:ncol(reverseECF)) {
    maxx <- ifelse(test = (sum(reverseECF[, col]) != 0), 
                   col, maxx)
    maxy <- max(which(reverseECF[, col] != 0), maxy)
  }
  df <- data.frame(Internesting.Period = numeric())
  for (col in 1:maxy) df <- cbind(df, A = numeric())
  colnames(df) <- c("Internesting.Period", paste0("Clutch.", 
                                                  0:(maxy - 1)))
  df0 <- df
  for (rowi in IP) {
    df0[1, ] <- 0
    df0[1, ] <- unname(c(rowi, reverseECF[1:maxy, rowi + 
                                            1]))
    df <- rbind(df, df0)
  }
  return(df)
}
