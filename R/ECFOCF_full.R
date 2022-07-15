#' ECFOCF_full calculate a table of probabilities of ECF and OCF.
#' @title Calculate a table of probabilities of ECF and OCF.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Return a matrix of class TableECFOCF.\cr
#' @param mu The average of lognormal for clutch frequency.
#' @param sd The sd parameter of lognormal for clutch frequency.
#' @param p The capture probability for an individual nesting event. As a logit
#' @param a The common capture probability. As a probability
#' @param OTN The relative probability of categories
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
#' If mu_season and sd_season are not NA, the model returns a 3D-table OCFECF.\cr
#' @family Model of Clutch Frequency
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' 
#' modelECFOCF <- ECFOCF_full(mu=c(mu1=5.58013243236187), 
#'                     sd=c(sd1=1.225581130238), 
#'                     mu_season=c(mu_season1=12), 
#'                     sd_season=c(sd_season1=2), 
#'                     p=c(p1=logit(0.7954041)), 
#'                     a=c(a1=1), 
#'                     MaxNests=15, 
#'                     MeanDaysBetween2Nests=9.8, 
#'                     length_season=floor(365/9.8)+1)
#' plot(modelECFOCF, period=12, max.scale=0.02)
#' modelECFOCF <- ECFOCF_full(mu=c(mu1=5.58013243236187), 
#'                     sd=c(sd1=1.225581130238), 
#'                     a=c(a1=1), 
#'                     p=c(p1=invlogit(1.3578137414575)), 
#'                     MaxNests=15)
#' plot(modelECFOCF)
#' }
#' @export

# Calcul table ECF OCF ####

ECFOCF_full <- function(mu, sd = NA, p, a=NULL, MaxNests=15, 
                     mu_season=NA, sd_season=NA, 
                     OTN=c(OTN1=1), 
                     MeanDaysBetween2Nests=9.8, 
                     length_season=NA, 
                     parallel=TRUE) {
  # mu <- NULL; sd <- NULL; p <- NULL; a <- NULL; MaxNests <- 15 
  # OTN=c(OTN1=1)
  # mu_season  <- NA; sd_season <- NA 
  # MeanDaysBetween2Nests <- 9.8
  # length_season  <- NA
  # parallel <- TRUE
  
  mln <- length(OTN)
  
  if (is.na(length_season)) {
    d3 <- 1
  } else {
    d3 <- length_season+MaxNests+1
  }
  
  dim_OCFECF <- c(MaxNests+1, MaxNests+1, d3)
  
  OCFECF <- array(data = 0, dim=dim_OCFECF, 
                  dimnames = list(paste0("OCF", 0:MaxNests), 
                                  paste0("ECF", 0:MaxNests), 
                                  paste0("time", 1:d3)))
  
  if (is.null(a)) {
    a <- rep(1, times=mln)
    names(a) <- paste0("a", 1:mln)
  }
  if (is.null(names(a))) names(a) <- paste0("a", 1:length(a))
  names(a)[grepl("a$", names(a))] <- gsub("a", "a1", names(a)[grepl("a$", names(a))])
  
  if (is.null(names(OTN))) names(OTN) <- paste0("OTN", 1:length(OTN))
  names(OTN)[grepl("OTN$", names(OTN))] <- gsub("OTN", "OTN1", names(OTN)[grepl("OTN$", names(OTN))])
  
  if (is.null(names(mu))) names(mu) <- paste0("mu", 1:length(mu))
  names(mu)[grepl("mu$", names(mu))] <- gsub("mu", "mu1", names(mu)[grepl("mu$", names(mu))])
  names(mu)[grepl("mu\\.", names(mu))] <- gsub("mu", "mu1", names(mu)[grepl("mu\\.", names(mu))])
  
  if (is.null(names(sd))) names(sd) <- paste0("sd", 1:length(sd))
  names(sd)[grepl("sd$", names(sd))] <- gsub("sd", "sd1", names(sd)[grepl("sd$", names(sd))])
  
  if (is.null(names(mu_season))) names(mu_season) <- paste0("mu_season", 1:length(mu_season))
  names(mu_season)[grepl("mu_season$", names(mu_season))] <- gsub("mu_season", "mu_season1", names(mu_season)[grepl("mu_season$", names(mu_season))])
  
  if (is.null(names(sd_season))) names(sd_season) <- paste0("sd_season", 1:length(sd_season))
  names(sd_season)[grepl("sd_season$", names(sd_season))] <- gsub("sd_season", "sd_season1", names(sd_season)[grepl("sd_season$", names(sd_season))])
  
  
  
  for (j in 1:mln) {
    p_ec <- p
    # je dois aussi prendre des p - 27-03-2021
    names(p_ec) <- gsub("p\\.", paste0("p", as.character(j), "."), names(p_ec))
    if (!all(identical(which(names(p_ec) == "p"), integer(0)))) {
      names(p_ec)[which(names(p_ec) == "p")] <- paste0("p", as.character(j))
    }
    
    pcommon <- p_ec[(names(p_ec)=="p") | (names(p_ec) == paste0("p", as.character(j)))][1]
    # if (is.na(pcommon)) pcommon <- NA
    if (dim_OCFECF[3]>1) {
      # p_period <- structure(rep(pcommon, dim(data)[3]-MaxNests),
      #                     .Names=paste0("p", as.character(j), ".",
      #                                   formatC(1:(dim(data)[3]-MaxNests), width=2, flag="0")))
      p_period <- structure(rep(pcommon, dim_OCFECF[3]),
                            .Names=paste0("p", as.character(j), ".",
                                          formatC(1:(dim_OCFECF[3]), width=2, flag="0")))
      
    } else {
      p_period <- structure(pcommon,
                            .Names=paste0("p", as.character(j)))
    }
    
    m1 <- match(names(p_ec), names(p_period))
    # if (all(is.na(m1))) m1 <- match(names(p_ec), paste0("p", as.character(j), ".",
    #              formatC(1:(dim(data)[3]-MaxNests), width=2, flag="0")))
    if (all(is.na(m1))) m1 <- match(names(p_ec), paste0("p", as.character(j), ".",
                                                        formatC(1:(dim_OCFECF[3]), width=2, flag="0")))
    p_period[m1[!is.na(m1)]] <- p_ec[!is.na(m1)]
    p_period <- 1/(1+exp(-p_period))
    
    p_period[m1[!is.na(m1)]] <- p_period[m1[!is.na(m1)]] * a[paste0("a", j)]
    
    nm <- paste0("p", as.character(j))
    
    OCFECF_ec <- ECFOCF_f(mu=mu[grepl(paste0("mu", j), names(mu))], 
                          sd=sd[paste0("sd", j)], 
                          p=p_period[substr(names(p_period), 1, nchar(nm))==nm],
                          mu_season = mu_season[paste0("mu_season", j)], 
                          sd_season = sd_season[paste0("sd_season", j)], 
                          MaxNests=MaxNests, 
                          length_season=length_season, 
                          parallel=parallel) * OTN[paste0("OTN", j)]
    
    
    OCFECF <- OCFECF+ OCFECF_ec
  }
    
    
    
    OCFECF <- addS3Class(OCFECF, "TableECFOCF")
    return(OCFECF)
}
