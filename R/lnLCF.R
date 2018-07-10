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
#' par <- c(par, fp[attributes(ECFOCF_2002)$table["begin"]:attributes(ECFOCF_2002)$table["end"]])
#' fixed.parameters <- c(p=-Inf)
#' 
#' lnLCF(x=par, data=ECFOCF_2002, fixed.parameters=fixed.parameters)
#' }
#' 
#' @export

# Log likelihood ECF OCF ####

lnLCF <- function(x, data, fixed.parameters=NULL, parallel=TRUE, verbose=FALSE) {
  
  x <- c(x, fixed.parameters)
  
  if (verbose) d(x)
  
  # dans ml j'ai le nombre max de catégories
  # La partie entière c'est la catégorie
  # Donc je peux créer un mu1.1
  ml <- floor(as.numeric(gsub("[a-zA-z]+", "", names(x))))
  if ((length(ml) == 1) | all(is.na(ml)) | (max(c(0, ml), na.rm=TRUE)==0)) {
    mln <- 1
  } else {
    mln <- max(ml, na.rm=TRUE)
  }
  
  MaxNests <- max(dim(data)[c(1, 2)])-1
  
  mu <- x[(substr(names(x), 1, 2)=="mu") & (substr(names(x), 1, 3) != "mu_")]
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
  
  sd <- x[(substr(names(x), 1, 2)=="sd") & (substr(names(x), 1, 3) != "sd_")]
  if (identical(sd, structure(numeric(0), .Names = character(0)))) sd <- c(sd=NA)
  if (length(sd)>1) sd <- sd[order(as.numeric(gsub("sd([0-9]+)", "\\1", names(sd))))]
  sd <- structure(c(sd, rep(sd[length(sd)], mln-length(sd))), .Names=paste0("sd", 1:mln))
  
  mu_season <- x[substr(names(x), 1, 9)=="mu_season"]
  if (length(mu_season)>1) mu_season <- mu_season[order(as.numeric(gsub("mu_season([0-9]+)", "\\1", names(mu_season))))]
  if (!identical(mu_season, structure(numeric(0), .Names = character(0)))) {
    mu_season <- c(mu_season, rep(mu_season[length(mu_season)], mln-length(mu_season)))
    names(mu_season) <- paste0("mu_season", 1:mln)
  } else {
    mu_season <- NA
  }
  
  sd_season <- x[substr(names(x), 1, 9)=="sd_season"]
  if (length(sd_season)>1) sd_season <- sd_season[order(as.numeric(gsub("sd_season([0-9]+)", "\\1", names(sd_season))))]
  if (!identical(sd_season, structure(numeric(0), .Names = character(0)))) {
    sd_season <- c(sd_season, rep(sd_season[length(sd_season)], mln-length(sd_season)))
    names(sd_season) <- paste0("sd_season", 1:mln)
  } else {
    sd_season <- NA
  }
  
  a <- x[substr(names(x), 1, 1)=="a"]
  if (identical(a, structure(numeric(0), .Names = character(0)))) {
    a <- structure(rep(Inf, mln), .Names=paste0("a", 1:mln))
  }
  if (any(names(a) == "a")) names(a[names(a) == "a"]) <- "a1"
  a_int <- structure(rep(Inf, mln), .Names=paste0("a", 1:mln))
  a_int[names(a)] <- a
  a <- 1/(1 + exp(-a_int))
  
  if (mln>1) {
    OTN <- abs(x[substr(names(x), 1, 3)=="OTN"])
    if (length(OTN)>1) OTN <- OTN[order(as.numeric(gsub("OTN([0-9]+)", "\\1", names(OTN))))]
    OTN <- c(OTN, 1)
    OTN <- c(OTN, rep(OTN[length(OTN)], mln-length(OTN)))
    names(OTN) <- paste0("OTN", 1:mln)
    OTN <- OTN/sum(OTN)
  } else {
    OTN <- c(OTN1=1)
  }
  
  p <- x[substr(names(x), 1, 1)=="p"]
  
  OCFECF <- data
  OCFECF[] <- 0
  
  for (j in 1:mln) {
    
    pcommon <- p[(names(p)=="p") | (names(p) == paste0("p", as.character(j)))][1]
    # if (is.na(pcommon)) pcommon <- NA
    if (dim(data)[3]>1) {
      p_period <- structure(rep(pcommon, dim(data)[3]-MaxNests),
                          .Names=paste0("p", as.character(j), ".",
                                        formatC(1:(dim(data)[3]-MaxNests), width=2, flag="0")))
    } else {
      p_period <- structure(pcommon,
                            .Names=paste0("p", as.character(j)))
    }

    m1 <- match(names(p), names(p_period))
    if (all(is.na(m1))) m1 <- match(names(p), paste0("p.",
                 formatC(1:(dim(data)[3]-MaxNests), width=2, flag="0")))

    p_period[m1[!is.na(m1)]] <- p[!is.na(m1)]
    p_period <- 1/(1+exp(-p_period))
    
    p_period[m1[!is.na(m1)]] <- p_period[m1[!is.na(m1)]] * a[paste0("a", j)]

    nm <- paste0("p", as.character(j))
    
    OCFECF <- OCFECF+ ECFOCF_f(mu=mu[grepl(paste0("mu", j), names(mu))], 
                               sd=sd[paste0("sd", j)], 
                               p=p_period[substr(names(p_period), 1, nchar(nm))==nm],
                               mu_season = mu_season[paste0("mu_season", j)], 
                               sd_season = sd_season[paste0("sd_season", j)], 
                               MaxNests=MaxNests, 
                               length_season=dim(data)[3]-MaxNests, 
                               parallel=parallel) * OTN[paste0("OTN", j)]
  }
  
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

