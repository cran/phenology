#' Tagloss_LengthObs returns a list with the number of days for different kinds of individuals are seen.
#' @title Return a list with the number of days for different kinds of individuals are seen.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Return a list with the number of days for different kinds of individuals are seen.
#' @param data Set of indivuals
#' @param progressbar Is shown a progressbar?
#' @description Usefull to summarize data
#' @family Model of Tag-loss
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' data_f_21 <- Tagloss_format(outLR, model="21")
#' Tagloss_LengthObs(data_f_21)
#' }
#' @export

Tagloss_LengthObs <- function(data, progressbar=TRUE) {
  
  lv <- levels(data$ID)
  
  if (progressbar) {
    if (is.element('progress', installed.packages()[,1])) {
      # library("progress")
      pb <- getFromNamespace("progress_bar", ns="progress")$new(
        format = "  completion [:bar] :percent eta: :eta",
        total = length(lv), clear = FALSE)
      libp <- TRUE
    } else {
      libp <- FALSE
      pb <- txtProgressBar(min=0, max=length(lv), style=3)
    }
  }
  
  L <- NULL
  R <- NULL
  cpt2 <- NULL
  cpt1 <- NULL
  cpt1R <- NULL
  cpt1L <- NULL
  
  cptlv <- 0
  
  for (IDx in lv) {
    
    if (progressbar) {
      if (libp) {
        pb$tick() 
        } else {
        cptlv <- cptlv + 1
        setTxtProgressBar(pb, cptlv)
      }
    }
    
    # data_ID <- subset(x = data, subset=(ID==IDx))
    data_ID <-data[data$ID==IDx, ]
    
    nL <- max(which(data_ID$L != ""))
    L <- c(L, data_ID$Date[nL] - data_ID$Date[1] + 1)
    nR <- max(which(data_ID$R != ""))
    R <- c(R, data_ID$Date[nR] - data_ID$Date[1] + 1)
    
    n2 <- which((data_ID$L != "") & (data_ID$R != ""))
    n1 <- which(((data_ID$L != "") & (data_ID$R == "")) | ((data_ID$L == "") & (data_ID$R != "")))
    n1L <- which(((data_ID$L != "") & (data_ID$R == "")))
    n1R <- which(((data_ID$L == "") & (data_ID$R != "")))
    
    if (identical(n2, integer(0))) {
      cpt2 <- c(cpt2, NA)
    } else {
      n2 <- max(n2)
      cpt2 <- c(cpt2, data_ID$Date[n2] - data_ID$Date[1] + 1)
    }
    if (identical(n1, integer(0))) {
      cpt1 <- c(cpt1, NA)
    } else {
      n1 <- max(n1)
      cpt1 <- c(cpt1, data_ID$Date[n1] - data_ID$Date[1] + 1)
    }
    if (identical(n1L, integer(0))) {
      cpt1L <- c(cpt1L, NA)
    } else {
      n1L <- max(n1L)
      cpt1L <- c(cpt1L, data_ID$Date[n1L] - data_ID$Date[1] + 1)
    }
    if (identical(n1R, integer(0))) {
      cpt1R <- c(cpt1R, NA)
    } else {
      n1R <- max(n1R)
      cpt1R <- c(cpt1R, data_ID$Date[n1R] - data_ID$Date[1] + 1)
    }
  }
  
  denR <- hist(log(R, base=10), plot=FALSE, breaks = seq(from=0, to=4, by=0.5))
  denL <- hist(log(L, base=10), plot=FALSE, breaks = seq(from=0, to=4, by=0.5))
  den2 <- hist(log(na.omit(cpt2), base=10), plot=FALSE, breaks = seq(from=0, to=4, by=0.5))
  den1 <- hist(log(na.omit(cpt1), base=10), plot=FALSE, breaks = seq(from=0, to=4, by=0.5))
  den1R <- hist(log(na.omit(cpt1R), base=10), plot=FALSE, breaks = seq(from=0, to=4, by=0.5))
  den1L <- hist(log(na.omit(cpt1L), base=10), plot=FALSE, breaks = seq(from=0, to=4, by=0.5))
  
  return(list(R=denR, L=denL, n2=den2, n1=den1, n1R=den1R, n1L=den1L))
}
