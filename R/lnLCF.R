#' lnLCF calculate the -log likelihood of data within a model.
#' @title Calculate the -log likelihood of data within a model.
#' @author Marc Girondot
#' @return Return the -log likelihood of data within a model.\cr
#' @param x A named vector of parameters (mu, sd, mu_season, sd_season, a, p and OTN).
#' @param data CMR database formated using TableECFOCF().
#' @param fixed.parameters Parameters that are fixed.
#' @param parallel If TRUE, parallel computing in ECFOCF_f is used.
#' @param verbose if TRUE, show the parameters.
#' @description Calculate the -log likelihood of data within a model.\cr
#' @family Model of Clutch Frequency
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' ECFOCF_2002 <- TableECFOCF(MarineTurtles_2002)
#' lnLCF(x=c(mu=4.71768454279272, 
#'                          sd=1.075711951667, 
#'                          p=-1.79746277312909), 
#'                  data=ECFOCF_2002)
#' 
#' 
#' ECFOCF_2002 <- TableECFOCF(MarineTurtles_2002, date0=as.Date("2002-01-01"))
#' fp <- rep(0, dim(ECFOCF_2002)[3])
#' names(fp) <- paste0("p.", formatC(1:(dim(ECFOCF_2002)[3]), width=2, flag="0"))
#' par <- c(mu1 = 0.6404831115214353, 
#'          sd1 = 0.69362774786433479, 
#'          mu2 = 5.6404831115214353, 
#'          sd2 = 5.69362774786433479, 
#'          mu_season = 12.6404831115214353, 
#'          sd_season = 1.69362774786433479, 
#'          OTN=1)
#' par <- c(par, fp[attributes(ECFOCF_2002)$table["begin"]:attributes(ECFOCF_2002)$table["final"]])
#' fixed.parameters <- c(p=-Inf)
#' 
#' lnLCF(x=par, data=ECFOCF_2002, fixed.parameters=fixed.parameters)
#' }
#' 
#' @export

# Log likelihood ECF OCF ####

lnLCF <- function(x, data, fixed.parameters=NULL, parallel=TRUE, verbose=FALSE) {
  
 #  x <- NULL; data <- NULL; fixed.parameters <- NULL; parallel <- TRUE; verbose <- TRUE
  
  xx <- c(x, fixed.parameters)
  
  if (verbose) d(xx)
  
  # dans ml j'ai le nombre max de catégories
  # La partie entière c'est la catégorie
  # Donc je peux créer un mu1.1
  ml <- suppressWarnings(floor(as.numeric(gsub("[a-zA-Z_]+", "", names(xx)))))
  if ((length(ml) == 1) | all(is.na(ml)) | (max(c(0, ml), na.rm=TRUE)==0)) {
    mln <- 1
  } else {
    mln <- max(ml, na.rm=TRUE)
  }
  
  MaxNests <- max(dim(data)[c(1, 2)])-1
  
  mu <- xx[(substr(names(xx), 1, 2)=="mu") & (substr(names(xx), 1, 3) != "mu_")]
  # if (length(mu)>1) mu <- mu[order(as.numeric(gsub("mu([0-9\\.]+)", "\\1", names(mu))))]
  # Ca ne va pas avec .
  # j'ai un mu. ou un mu
  # mu_ref <- NULL
  if (mln > 1) {
    if (any(names(mu)=="mu")) {
      mu_ref <- mu[names(mu)=="mu"]
      mu_ec <- NULL
      for (i in 1:mln) {
        if (all(!grepl(paste0("mu", i), names(mu)))) {
          mu_ec <- c(mu_ec, structure(unname(mu_ref), .Names=paste0("mu", i)))
        } else {
          mu_ec <- c(mu_ec, mu[grepl(paste0("mu", i), names(mu))])
        }
      }
      mu <- mu_ec
    }
    if (any(grepl("mu\\.+", names(mu)))) {
      mu_ref <- mu[grepl("mu\\.+", names(mu))]
      
      mu_ec <- NULL
      for (i in 1:mln) {
        if (all(!grepl(paste0("mu", i), names(mu)))) {
          mu_ec <- c(mu_ec, structure(unname(mu_ref), .Names=gsub("mu(\\.[0-9]+)", paste0("mu", i, "\\1"), names(mu_ref))))
        } else {
          mu_ec <- c(mu_ec, mu[grepl(paste0("mu", i), names(mu))])
        }
      }
      mu <- mu_ec
    }
  } else {
    names(mu) <- gsub("mu[0-9]*(\\.*[0-9]*)", "mu1\\1", names(mu))
  }
  
  sd <- xx[(substr(names(xx), 1, 2)=="sd") & (substr(names(xx), 1, 3) != "sd_")]
  if (identical(sd, structure(numeric(0), .Names = character(0)))) sd <- c(sd=NA)
  if (length(sd)>1) sd <- sd[order(as.numeric(gsub("sd([0-9]+)", "\\1", names(sd))))]
  sd <- structure(c(sd, rep(sd[length(sd)], mln-length(sd))), .Names=paste0("sd", 1:mln))
  
  mu_season <- xx[substr(names(xx), 1, 9)=="mu_season"]
  if (length(mu_season)>1) mu_season <- mu_season[order(as.numeric(gsub("mu_season([0-9]+)", "\\1", names(mu_season))))]
  if (!identical(mu_season, structure(numeric(0), .Names = character(0)))) {
    mu_season <- c(mu_season, rep(mu_season[length(mu_season)], mln-length(mu_season)))
    names(mu_season) <- paste0("mu_season", 1:mln)
  } else {
    mu_season <- NA
  }
  
  sd_season <- xx[substr(names(xx), 1, 9)=="sd_season"]
  if (length(sd_season)>1) sd_season <- sd_season[order(as.numeric(gsub("sd_season([0-9]+)", "\\1", names(sd_season))))]
  if (!identical(sd_season, structure(numeric(0), .Names = character(0)))) {
    sd_season <- c(sd_season, rep(sd_season[length(sd_season)], mln-length(sd_season)))
    names(sd_season) <- paste0("sd_season", 1:mln)
  } else {
    sd_season <- NA
  }
  
  a <- xx[substr(names(xx), 1, 1)=="a"]
  if (identical(a, structure(numeric(0), .Names = character(0)))) {
    a <- structure(rep(Inf, mln), .Names=paste0("a", 1:mln))
  }
  if (any(names(a) == "a")) names(a[names(a) == "a"]) <- "a1"
  a_int <- structure(rep(Inf, mln), .Names=paste0("a", 1:mln))
  a_int[names(a)] <- a
  a <- 1/(1 + exp(-a_int))
  
  if (mln>1) {
    OTN <- abs(xx[substr(names(xx), 1, 3)=="OTN"])
    if (length(OTN)>1) OTN <- OTN[order(as.numeric(gsub("OTN([0-9]+)", "\\1", names(OTN))))]
    OTN <- c(OTN, 1)
    OTN <- c(OTN, rep(OTN[length(OTN)], mln-length(OTN)))
    names(OTN) <- paste0("OTN", 1:mln)
    OTN <- OTN/sum(OTN)
  } else {
    OTN <- c(OTN1=1)
  }
  
  p <- xx[substr(names(xx), 1, 1)=="p"]
  # p <- 1/(1 + exp(-p))
  
  if (any((is.na(mu_season)) | (is.na((sd_season))))) {
    length_mean_ec <- NA
  } else {
    length_mean_ec  <-  attributes(data)$characteristics["length_season"]
  }
  
  
  OCFECF <- ECFOCF_full(mu=mu, 
                        sd=sd, 
                        p=p, 
                        a=a, 
                        MaxNests=MaxNests, 
                        mu_season=mu_season, 
                        sd_season=sd_season, 
                        OTN=OTN, 
                        MeanDaysBetween2Nests=attributes(data)$characteristics["MeanDaysBetween2Nest"], 
                        length_season=length_mean_ec, 
                        parallel=parallel)
  
  if (any(is.na(OCFECF[1, 1, ]))) {
    # ss1 <<- x
    return(+Inf)
  } else {
    # OCFECF <<- OCFECF
    # OCFECF[OCFECF >= 1-(1E-100)] <- 1-(1E-100)
    # OCFECF[OCFECF <= (1E-100)] <- (1E-100)
    
    OCF0 <- sum(OCFECF[1, 1, ])
    
    
    
    OCFECF <- OCFECF[-1, -1, ]/(1-OCF0)
    # OCFECF[OCFECF >= 1-(1E-100)] <- 1-(1E-100)
    # OCFECF[OCFECF <= (1E-100)] <- (1E-100)
    
    OCFECF <- log(OCFECF)
    OCFECF <- OCFECF*data[-1, -1, ]
    LnL <- -sum(OCFECF, na.rm = TRUE) 
    # ss2 <<- x
    return(LnL)
  }
}

