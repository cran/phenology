#' ECFOCF_f calculate a table of probabilities of ECF and OCF.
#' @title Calculate a table of probabilities of ECF and OCF.
#' @author Marc Girondot
#' @return Return a matrix of class TableECFOCF.\cr
#' @param mu The average of lognormal for clutch frequency.
#' @param sd The sd parameter of lognormal for clutch frequency.
#' @param p The capture probability for an individual nesting event.
#' @param MaxNests Maximum number of nests by a female.
#' @param MeanDaysBetween2Nests Average number of days between two nests.
#' @param mu_season The average of ordinal day for beginning of nesting season.
#' @param sd_season The sd parameter of lognormal for ordinal day for beginning of nesting season.
#' @param length_season The total length of season based on groups of interclutch intervals.
#' @param parallel If TRUE parallel computing is used.
#' @description This function calculates a table of probabilities of ECF and OCF.\cr
#' If p is lower or higher than 1E-100 or 1-1E-100, it is changed to 1E-100 and 1-(1E-100) respectively.\cr
#' Names for p vector elements should be p, or px (with x=1:categories), or px.period.\cr
#' If mu_season and sd_season are equal to NA, the model is not temporalized.\cr
#' If mu_season and sd_season are not NA, the model is temporalized.\cr
#' @family Model of Clutch Frequency
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' modelECFOCF <- ECFOCF_f(mu=5.58013243236187, 
#'                     sd=1.225581130238, 
#'                     p=invlogit(1.3578137414575), 
#'                     MaxNests=15)
#' plot(modelECFOCF)
#' modelECFOCF <- ECFOCF_f(mu=5.58013243236187, 
#'                     sd=1.225581130238, 
#'                     mu_season=12, 
#'                     sd_season=2, 
#'                     p=c(p1=invlogit(1.3578137414575)), 
#'                     MaxNests=15, 
#'                     MeanDaysBetween2Nests=9.8, 
#'                     length_season=floor(365/9.8)+1
#'                     )
#' plot(modelECFOCF, period=2)
#' }
#' @export

# Calcul table ECF OCF ####

ECFOCF_f <- function(mu, sd, p, MaxNests=15, 
                     mu_season=NA, sd_season=NA, 
                     MeanDaysBetween2Nests=9.8, 
                     length_season=floor(365/MeanDaysBetween2Nests)+1, 
                     parallel=TRUE) {
  
  
  if (is.na(mu_season) | is.na(sd_season)) {
    x <- 1:MaxNests
    p <- ifelse(p < 1E-100, 1E-100, p)
    p <- ifelse(p > 1 - 1E-100, 1 - 1E-100, p)
    
    p <- p[1]
    
    OCFECF <- array(data = 0, dim=c(MaxNests+1, MaxNests+1, 1), 
                    dimnames = list(paste0("OCF", 0:(MaxNests)), paste0("ECF", 0:(MaxNests)), 
                                    "time1"))
 
  y <- dlnorm(x, meanlog=log(abs(mu)), sdlog=abs(sd))
  y <- y / sum(y)
  for (x_ec in x) {
    # Dans x_ec, j'ai le CF
    # Dans OCF_ec, j'ai de 0 à CF, la probabilité d'observer rang+1 pontes
    diverses <- expand.grid(rep(list(0:1), x_ec))
    OCF_ec <- apply(diverses, MARGIN = 1, FUN = sum)
    ECF_ec <- apply(diverses, MARGIN = 1, FUN = function(ec) {
      if (sum(ec) != 0) {
        max(which(ec == 1))-min(which(ec == 1))+1
      } else {0}
    })
    
    p_ec <- p^OCF_ec*(1-p)^(x_ec-OCF_ec)
    for (i in seq_along(ECF_ec)) {
      OCFECF[OCF_ec[i]+1, ECF_ec[i]+1, 1] <- (p_ec[i]*y[x_ec]) + OCFECF[OCF_ec[i]+1, ECF_ec[i]+1, 1]
    }
  }

  } else {
    # Je suis dans un modèle 3D
    x <- 1:MaxNests
    nm <- names(p)[1]
    if (grepl("\\.", names(p)[1])) nm <- gsub("(p[0-9])+\\.[0-9]+$", "\\1", names(p)[1])
    nm <- gsub("p([0-9])+$", "\\1", nm)
    nm <- ifelse(nm=="p", 1, as.numeric(nm))
    
    # je crée une chaîne avec toute les prob nommées
    commonprob <- ifelse(is.na(p[paste0("p", as.character(nm))]), 0, p[paste0("p", as.character(nm))])
    prob_ec <- structure(rep(commonprob, length_season+MaxNests), 
                         .Names=paste0("p", as.character(nm), ".", formatC(1:(length_season+MaxNests), width=2, flag="0")))
    
    
    # Oui je met les valeurs à leur place
    prob_ec[names(p)[grepl("\\.", names(p))]] <- p[names(p)[grepl("\\.", names(p))]]

    p <- prob_ec
    
    # p <- ifelse(p < 1E-9, 1E-9, p)
    # p <- ifelse(p > 1 - 1E-9, 1 - 1E-9, p)
    
    # length_season+MaxNests pas length_season+MaxNests+1
    OCFECF <- array(data = 0, dim=c(MaxNests+1, MaxNests+1, length_season+MaxNests), 
                    dimnames = list(paste0("OCF", 0:(MaxNests)), paste0("ECF", 0:(MaxNests)), 
                                    # Là aussi
                                    paste0("time", 1:(length_season+MaxNests))
                                    )
                    )
    
    time <- dlnorm(1:length_season, meanlog=log(abs(mu_season)), 
                   sdlog=abs(sd_season))
    time <- time / sum(time)
    
    # j'ai length_season valeurs et pour chacune la probabilité
    
    if ((.Platform$OS.type == "unix") & (parallel)) {
      cores <- detectCores()
    } else {
      cores <- 1
    }
    
    tlist <- list()
    if (cores > 1) {
      for (i in 1:(cores-1)) tlist <- c(tlist, list(seq(from=1, to=floor(length_season/cores), by=1)+((i-1)*floor(length_season/cores))))
      tlist <- c(tlist, list(seq(from=rev(tlist[[cores-1]])[1]+1, to=length_season, by=1)))
    } else {
      tlist <- list(seq(from=1, to=length_season, by=1))
    }
    
    y <- dlnorm(x, meanlog=log(abs(mu)), sdlog=abs(sd))
    y <- y / sum(y)
    
    # y[1] correspond à la valeur x=1 donc O car x-1
    
    diverses_list <- list()
    first_obs_list <- list()
    OCF_ec_list <- list()
    ECF_ec_list <- list()
    
    for (x_ec in x) {
      # Dans x_ec, j'ai le CF
      # Dans OCF_ec, j'ai de 0 à CF, la probabilité d'observer rang+1 pontes
      diverses <- expand.grid(rep(list(0:1), x_ec))
      diverses_list <- c(diverses_list, list(diverses))
      
      first_obs <- apply(diverses, MARGIN = 1, FUN = function(ec) {
        if (sum(ec) != 0) {
          min(which(ec == 1))-1
        } else {0}
      })
      first_obs_list <- c(first_obs_list, list(first_obs))
      
      
      OCF_ec <- apply(diverses, MARGIN = 1, FUN = sum)
      OCF_ec_list <- c(OCF_ec_list, list(OCF_ec))
      
      ECF_ec <- apply(diverses, MARGIN = 1, FUN = function(ec) {
        if (sum(ec) != 0) {
          max(which(ec == 1))-min(which(ec == 1))+1
        } else {0}
      })
      ECF_ec_list <- c(ECF_ec_list, list(ECF_ec))
    }
    
   
    OCFECF_L <- mclapply(X=tlist, mc.cores = cores, 
                         FUN=function(t_int) {
      OCFECF_int <- OCFECF
      # print(t)
      
      
      for (t in t_int) {
      
      # dans la la distribution des CF
      # ils commencent tous le jour time
      # mais certains peuvent ne pas avoir été vu
      
      for (x_ec in x) {
        # Dans x_ec, j'ai le CF
        # Dans OCF_ec, j'ai de 0 à CF, la probabilité d'observer rang+1 pontes
        diverses <- diverses_list[[x_ec]]
        
        ptcf <- y[x_ec] * time[t]
        
        first_obs <- first_obs_list[[x_ec]]
        OCF_ec <- OCF_ec_list[[x_ec]]
        ECF_ec <- ECF_ec_list[[x_ec]]
        
        saec <- t:(t+x_ec-1)
        p_ec <- apply(diverses, MARGIN = 1, FUN = function(ec) {
          prod(p[saec]^ec*(1-p[saec])^(1-ec))
        })
        p_ec <- p_ec * ptcf
        
        for (i in seq_along(ECF_ec)) {
          OCFECF_int[OCF_ec[i]+1, ECF_ec[i]+1, t+first_obs[i]] <- 
            OCFECF_int[OCF_ec[i]+1, ECF_ec[i]+1, t+first_obs[i]] + p_ec[i]
        }
      }
      }
      return(OCFECF_int)
    })
    

    for (i in seq_along(OCFECF_L)) OCFECF <- OCFECF + OCFECF_L[[i]]
    
    # fin de calcul de l'array
  }
    
    
    # OCFECF <- OCFECF / (1-sum(OCFECF[1, 1, ]))
    class(OCFECF) <- "TableECFOCF"
    return(OCFECF)
}
