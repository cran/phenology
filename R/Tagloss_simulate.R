#' Tagloss_simulate returns a list with the number of days different kinds of individuals are seen.
#' @title Return a list with the number of days different kinds of individuals are seen.
#' @author Marc Girondot
#' @return Return a list with the number of days different kinds of individuals are seen.
#' @param n Number of individuals to simulate
#' @param par Set of parameters
#' @param pobservation Probability of daily observation
#' @param dailysurvival Daily probability of survival
#' @param LengthObservation The log of number of days of observations is a random number between min and max
#' @param model Must be "12" or "LR"
#' @param model_before Transformation of parameters before to use Tagloss_model()
#' @param model_after Transformation of parameters after to use Tagloss_model()
#' @param progressbar Is a progressbar should be shown?
#' @description Generate data with known features.\cr
#'   model_before is applied to par parameter.\cr
#'   model_after is applied after par is separated in p1, p2, pL1, pL2, pR1 and pR2 parameters.\cr
#' pobservation can be a vector of daily probabilities to be captured. The last value is repeated 
#' if necessary.\cr
#' The maximum number of days of observation is exp(LengthObservation["max"]).\cr
#' If model="12" then par must have _1 and _2 parameters.\cr
#' if model="LR" then par must have _L2, _L1, _R2, R1 parameters.\cr
#' @family Model of Tag-loss
#' @examples
#' library(phenology)
#' \dontrun{
#' # Example
#' par <- structure(c(49.5658922243074, 808.136085362158, 106.283783786853, 
#' 5.22150592456511, 8.00608716525864, 8.32718202233396, 150.612916258503, 
#' 715.865805125223, 2242.06574225966, 119.212383120678, 10.1860735529433, 
#' 7.14231725937626), .Names = c("D1_2", "D2D1_2", "D3D2_2", "A_2", 
#' "B_2", "C_2", "D1_1", "D2D1_1", "D3D2_1", "A_1", "B_1", "C_1"))
#' cmr <- Tagloss_simulate(n=500, 
#'                         par=par, model="12")
#' cmr_f <- Tagloss_format(cmr, model="12")
#' }
#' @export


# @param daiysurvival Probability of daily survivorship

Tagloss_simulate <- function(n=500, par, pobservation=c(rep(0.05, 70), 0.01), 
                             LengthObservation=c(min=0, max=9), 
                             dailysurvival=0.999, 
                             model="12", 
                             model_before=NULL, model_after=NULL, 
                             progressbar=TRUE) {
  # model_before is applied to par parameter
  # model_after is applied after par is separated in pL1, pL2, pR1 and pR2 parameters
  
  # dailysurvival=dailysurvival
  
  Lastpobs <- rev(x = pobservation)[1]
  lp <- length(pobservation)
  
  out <- data.frame(ID=character(), L=character(), R=character(), Date=as.Date(as.character()), 
                    stringsAsFactors = FALSE)
  if (progressbar) {
    if (is.element('progress', installed.packages()[,1])) {
      # library("progress")
      pb <- getFromNamespace("progress_bar", ns="progress")$new(
        format = "  completion [:bar] :percent eta: :eta",
        total = n, clear = FALSE)
      libp <- TRUE
    } else {
      libp <- FALSE
      pb <- txtProgressBar(min=0, max=n, style=3)
    }
  }
  
  if (!is.null(model_before)) eval(parse(text = model_before))
  pL1 <- par[grepl("_L1", names(par))]
  pL2 <- par[grepl("_L2", names(par))]
  pR1 <- par[grepl("_R1", names(par))]
  pR2 <- par[grepl("_R2", names(par))]
  p1 <- par[grepl("_1", names(par))]
  p2 <- par[grepl("_2", names(par))]
  if (!is.null(model_after)) eval(parse(text = model_after))
  
  if (!identical(unname(pL1), numeric(0))) pL1 <- Tagloss_model(1:floor(exp(LengthObservation["max"])), pL1)
  if (!identical(unname(pL2), numeric(0))) pL2 <- Tagloss_model(1:floor(exp(LengthObservation["max"])), pL2)
  if (!identical(unname(pR1), numeric(0))) pR1 <- Tagloss_model(1:floor(exp(LengthObservation["max"])), pR1)
  if (!identical(unname(pR2), numeric(0))) pR2 <- Tagloss_model(1:floor(exp(LengthObservation["max"])), pR2)
  if (!identical(unname(p1), numeric(0))) p1 <- Tagloss_model(1:floor(exp(LengthObservation["max"])), p1)
  if (!identical(unname(p2), numeric(0))) p2 <- Tagloss_model(1:floor(exp(LengthObservation["max"])), p2)
  
  if (model == "12") {
    pL1 <- p1
    pR1 <- p1
    pL2 <- p2
    pR2 <- p2
  }
  
  for (ID in 1:n) {
  if (progressbar) {
    if (libp) pb$tick() else setTxtProgressBar(pb, ID)
    }
    
    repeat {
      LObs <- floor(exp(runif(1, min=LengthObservation["min"], max=LengthObservation["max"])))
      
      # Je modifie LObs en fonction de dailysurvival
      LObsmax <- which(runif(LObs) > dailysurvival)
      if (!identical(LObsmax, integer(0))) LObs <- min(LObsmax)
    
    history <- data.frame(ID=as.character(ID), 
                          L=rep(paste0(as.character(ID), "_L"), LObs), 
                          R=rep(paste0(as.character(ID), "_R"), LObs), 
                          Date=seq(from=as.Date("2000-01-01"), length.out = LObs, by="1 day"), 
                          stringsAsFactors = FALSE)
    
    perteL2 <- which(runif(LObs) <= pL2[1:LObs])
    if (!identical(perteL2, integer(0))) perteL2 <- min(perteL2)
    perteR2 <- which(runif(LObs) <= pR2[1:LObs])
    if (!identical(perteR2, integer(0))) perteR2 <- min(perteR2)
    
    if ((!identical(perteL2, integer(0))) & (!identical(perteR2, integer(0)))) {
      # Il a perdu les deux; je prends la plus petite
      if (perteL2 <= perteR2) {
        # Il perd d'abord la L
        history[perteL2:LObs, "L"] <- ""
        perteR1 <- which(runif(LObs-perteL2+1)<=pR1[perteL2:LObs]) + perteL2 - 1
        if (!identical(perteR1, numeric(0))) history[min(perteR1):LObs, "R"] <- ""
      } else {
        history[perteR2:LObs, "R"] <- ""
        perteL1 <- which(runif(LObs-perteR2+1)<=pL1[perteR2:LObs]) + perteR2 - 1
        if (!identical(perteL1, numeric(0))) history[min(perteL1):LObs, "L"] <- ""
      }
    }
    if ((!identical(perteL2, integer(0))) & (identical(perteR2, integer(0)))) {
      # Il perd d'abord la L
      history[perteL2:LObs, "L"] <- ""
      perteR1 <- which(runif(LObs-perteL2+1)<=pR1[perteL2:LObs]) + perteL2 - 1
      if (!identical(perteR1, numeric(0))) history[min(perteR1):LObs, "R"] <- ""
    }
    if ((!identical(perteR2, integer(0))) & (identical(perteL2, integer(0)))) {
      # Il perd d'abord la R
      history[perteR2:LObs, "R"] <- ""
      perteL1 <- which(runif(LObs-perteR2+1)<=pL1[perteR2:LObs]) + perteR2 - 1
      if (!identical(perteL1, numeric(0))) history[min(perteL1):LObs, "L"] <- ""
    }
    
    # J'ai généré les histoires de perte des bagues
    # Maintenant je les vois quand ?
    history1 <- history[1, ]
    if (lp < LObs) pobs_encours <-  c(pobservation, rep(Lastpobs, LObs-lp)) else pobs_encours <- pobservation[1:LObs]
    history[which(runif(LObs) > pobs_encours), "ID"] <- NA
    history[1, ] <- history1
    history <- na.omit(history)
    
    if ((nrow(history) > 1) & (history[1, "L"] != "") & (history[1, "R"] != "")) break
    }

    
    out <- rbind(out, history)
    }
  out$ID <- as.factor(out$ID)
  return(out)
  
}