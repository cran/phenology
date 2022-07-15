#' generateCF generates set of data to test fitCF.
#' @title Generate a set of data to test Clutch Frequency for marine turtles.
#' @author Marc Girondot
#' @return Return a list with 4 elements: Category, CF, Beginning and Observations being a 
#' dataframe of individuals.
#' @param x Initial parameters to be used
#' @param MeanDaysBetween2Nests Number of days in average between two nests
#' @param date0 Initial date to generate data
#' @param n Number of individuals to model
#' @param verbose If TRUE, give information about each animal.
#' @description This function generates a dataframe to test \code{fitCF()}.\cr
#' This model is an enhanced version of the one published by Briane et al. (2007).\cr
#' Parameters are \code{mu} and \code{sd} being the parameters of a  
#' distribution used to model the clutch frequency.\cr
#' This distribution is used only as a guide but has not statistical meaning.\cr
#' The parameter \code{p} is the -logit probability that a female is seen 
#' on the beach for a particular nesting event. It includes both the probability 
#' that it is captured but also the probability that it uses that specific beach.\cr
#' Several categories of females can be included in the model using index after 
#' the name of the parameter, for example \code{mu1}, \code{sd1} and \code{mu2}, 
#' \code{sd2} indicates that two categories of females with different clutch 
#' frequencies distribution are present. Similarly \code{p1} and \code{p2} indicates 
#' that two categories of females with different capture probabilities are present.\cr
#' If more than one category is used, then it is necessary to include the 
#' parameter \code{OTN} to indicate the relative frequencies of each category. 
#' If two categories are used, one \code{OTN} parameter named \code{ONT1} must 
#' be included. The \code{OTN2} is forced to be 1. Then the relative frequency 
#' for category 1 is \code{OTN1/(OTN1+1)} and for category 2 is \code{1/(OTN1+1)}. 
#' Same logic must be applied for 3 and more categories with always the last one 
#' being fixed to 1.\cr
#' 
#' if p or a (logit of the capture probability) are equal to -Inf, 
#' the probability of capture is 0 and if they are equal to 
#' +Inf, the probability is 1.\cr
#' 
#' The value of p out of the period 
#' of nesting must be set to +Inf (capture probability=1)
#' to indicate that no turtle is nesting in this period.\cr
#' 
#' p must be set to -Inf (capture probability=0) to indicate that no
#' monitoring has been done during a specific period of the nesting season.\cr
#' 
#' The best way to indicate capture probability for 3D model (OCF, ECF, Period) 
#' is to indicate p.period common for all categories and a1, a2, etc for each category. 
#' The capture probability for category 1 will be p.period * a1, and for category 2 
#' will be p.period * a2, etc. \cr
#' 
#' In this case, the parameters p.period should be indicated in fitted parameters 
#' as well as a1, but a2 must be fixed to +Inf in fixed.parameters. Then the capture 
#' probability for category 2 will be p.period and for category 1 a1 * p.period.\cr
#' 
#' @family Model of Clutch Frequency
#' @seealso Briane J-P, Rivalan P, Girondot M (2007) The inverse problem applied 
#'             to the Observed Clutch Frequency of Leatherbacks from Yalimapo beach, 
#'             French Guiana. Chelonian Conservation and Biology 6:63-69
#' @seealso Fossette S, Kelle L, Girondot M, Goverse E, Hilterman ML, Verhage B, 
#'          Thoisy B, de, Georges J-Y (2008) The world's largest leatherback 
#'          rookeries: A review of conservation-oriented research in French 
#'          Guiana/Suriname and Gabon. Journal of Experimental Marine Biology 
#'          and Ecology 356:69-82
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' 
#' par <- c(mu = 2.4911638591178051, 
#'          sd = 0.96855483039640977, 
#'          mu_season = 13.836059118657793, 
#'          sd_season = 0.17440085345943984, 
#'          p.10 = 1.3348233607728222, 
#'          p.11 = 1.1960387774393837, 
#'          p.12 = 0.63025680979544774, 
#'          p.13 = 0.38648155002707452, 
#'          p.14 = 0.31547864054366048, 
#'          p.15 = 0.19720001827017075, 
#'          p.16 = 0.083199496372073328, 
#'          p.17 = 0.32969130595897905, 
#'          p.18 = 0.36582777525265819, 
#'          p.19 = 0.30301248314170637, 
#'          p.20 = 0.69993987591518514, 
#'          p.21 = 0.13642423871641118, 
#'          p.22 = -1.3949268190534629, 
#'          p=+Inf)
#' 
#' o_mu1p1season1 <- generateCF(x=par, n=1, verbose=TRUE)
#' o_mu1p1season1 <- generateCF(x=par, n=1000)
#' plot(o_mu1p1season1$CF)
#' hist(o_mu1p1season1$Beginning)
#' }
#' @export

generateCF <- function(x=c(mu=4, sd=1, 
                           p=+Inf, 
                           mu_season = 13.836059118657793, 
                           sd_season = 0.17440085345943984),
                       MeanDaysBetween2Nests=9.8, 
                       date0=as.Date("2020-01-01"), 
                       n=1, verbose = TRUE) {
  
  data <- data.frame(Index=character(), Date=as.Date(character()))
  OTN <- x[substr(names(x), 1, 3) == "OTN"]
  lnm <- max(c(1, floor(as.numeric(gsub("[A-Za-z_]", "", names(x))))), 
             na.rm=TRUE)
  if (identical(structure(numeric(0), .Names = character(0)), OTN)) {
    OTN <- c(OTN1=1)
  }
  
  if (verbose) print(paste0("I have detected ",as.character(lnm), " categor", ifelse(lnm>1, "ies.", "y.")))
  if (verbose) {
    print("The relative probabilities are:")
    print(paste0(names(OTN), "=", specify_decimal(unname(OTN), decimals=3)))
  }
  
  CF_tot <- NULL
  Beginning <- NULL
  cat <- NULL
  
  for (i in 1:n) {
    
    if (verbose) print(paste0("Individual ", as.character(i)))
    nest <- rep(0, 100)
    
    nm <- min(which(runif(1)<c(OTN, 1)))
    cat <- c(cat, nm)
    if (verbose) print(paste0("This is individual of category ",as.character(nm)))
    
 
    p_freq <- rep(0, 100)
    names(p_freq) <- paste0("p.", formatC(1:100, width = 2, flag = "0"))
    
    mu_season <- x[names(x) == "mu_season" | 
                     names(x) == paste0("mu_season", 
                                        formatC(nm, width=2, flag="0"))]
    sd_season <- x[names(x) == "sd_season" | 
                     names(x) == paste0("sd_season", 
                                        formatC(nm, width=2, flag="0"))]
    mu  <- x[names(x) == "mu" | 
               names(x) == paste0("mu", 
                                  formatC(nm, width=2, flag="0"))]

    sd <- x[names(x) == "sd" | 
                   names(x) == paste0("sd", 
                                      formatC(nm, width=2, flag="0"))]
    
    deb <- floor(rlnorm(1, meanlog = log(abs(mu_season)), 
           sdlog = abs(sd_season))+0.5)
    
    Beginning <- c(Beginning, deb)
    
    if (verbose) print(paste0("The period for first nest is ", deb))
    
    cf <- floor(rlnorm(1, meanlog=log(abs(mu)), sdlog=abs(sd))+0.5)
    
    CF_tot <- c(CF_tot, cf)
    
    if (verbose) print(paste0("The clutch frequency is ", cf))
    
    p <- x[substr(names(x), 1, 1) == "p" | 
              substr(names(x), 1, 1) == paste0("p", 
                                 formatC(nm, width=2, flag="0"))]
    names(p) <- gsub(paste0("p", formatC(nm, width=2, flag="0"), "([.$])"), 
                      paste0("p", "\1"), names(p))
    
    if (!is.na(p["p"])) p_freq[] <- unname(p["p"])
    
    
    p_index <- p[match(names(p_freq), names(p))]
    p_index <- p_index[!is.na(p_index)]
    
 
    p_freq[!is.na(match(names(p_freq), names(p)))] <- unname(p_index)
    p_freq <- invlogit(p_freq)

    for (pos in deb:(deb+cf-1)) {
      if (runif(1) < p_freq[pos]) nest[pos] <-1
    }
    
    if (any(nest != 0)) {
      
      data_int <- data.frame(Index=rep(as.character(i)), 
                             Date=date0 + (which(nest != 0)-1)*MeanDaysBetween2Nests)
      data <- rbind(data, data_int)
      if (verbose) print(nest)
    } else {
      if (verbose) print("No nest is observed.")
    }
  }
  
  data <- list(Category=cat, CF=CF_tot, 
               Beginning=Beginning, Observations=data)
  
  return(data)
}

