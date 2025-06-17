#' RI2BP calculate Breeding Proportion from Remigration Interval.
#' @title Calculate Breeding Proportion from Remigration Interval.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Return a vector with mean and se.\cr
#' @param proportion A vector of the proportion of different categories
#' @param RI The RI of individuals - See description
#' @param sampling Number of individuals for sampling
#' @param replicates Number of replicates to estimate SE
#' @description This function calculates breeding proportion (BP) from Remigration Interval (RI).\cr
#' RI can be a vector of RI, one per individual, or an aggregated value with mean and sd and 
#' several categories can exist. See examples.
#' @family Model of Remigration Interval
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' RI <- c(mean=3, sd=0)
#' RI2BP(RI=RI)
#' 
#' RI <- c(mean=2, sd=0.3)
#' RI2BP(RI=RI)
#' 
#' RI <- c(mean=4, sd=0)
#' RI2BP(RI=RI)
#' 
#' RI <- c(mean=c(2, 10), sd=c(0, 0))
#' proportion <- c(0.5, 0.5)
#' RI2BP(proportion=proportion, RI=RI)
#' 
#' c <- c(c1=0.1, c2=0.72, c3=0.126, c4=0.0378, c5=0.0162)
#' plot(1:5, c, xlab="RI", ylab="Proportion", las=1, bty="n", type="h")
#' RI2BP(proportion=c, RI=c(mean=1:5, sd=rep(0, 5)))
#' 
#' # To generate random RI with known mean and sd using 
#' # a truncated lognormal distribution (because 0 does not exist)
#' mean <- 2.5 - 1
#' sd <- 1
#' location <- log(mean^2 / sqrt(sd^2 + mean^2)) 
#' shape <- sqrt(log(1 + (sd^2 / mean^2)))   
#' RI <- round(rlnorm(100, meanlog = location, sdlog = shape))+1
#' mean(RI); sd(RI) # All is ok !
#' plot(table(RI), xlab="RI", ylab="Proportion", las=1, bty="n", type="h")
#' RI2BP(proportion=proportion, RI=RI)
#' 
#' }
#' @export

# Calcul table ECF OCF ####

RI2BP <- function(proportion=1, RI=c(mean=2, sd=0), sampling=10000, replicates=100) {
  
  if (any(grepl("mean", names(RI)))) {
    np <- length(RI) / 2
    if (np == 1) {
      proportion <- 1
      names(RI)[which(names(RI) == "mean")] <- "mean1"
      names(RI)[which(names(RI) == "sd")] <- "sd1"
      
    } else {
      if (length(proportion) == np - 1) {
        proportion <- c(proportion, 1 - sum(proportion))
      }
    }
    
    if (length(proportion) != np)
      stop("Length of proportion is not compatible with RI.")
    
    means <- RI[paste0("mean", as.character(1:np))]
    sds <- RI[paste0("sd", as.character(1:np))]
    
    means <- means - 1
    
    locations <- log(means^2 / sqrt(sds^2 + means^2))
    shapes <- sqrt(log(1 + (sds^2 / means^2)))
    
    zero <- (is.na(locations) & is.na(shapes))
    
    if (any(zero)) {
      locations[zero] <- 0
      shapes[zero] <- 0
    }
    # r <- round(rlnorm(n=1000000, location, shape)) + 1
    
    if (any(is.na(locations)) | any(is.na(shapes))) stop("Check the RI, something is wrong")
    
    rep_RIi <- NULL
    for (r in 1:replicates) {
      category <- sample(1:np, sampling, replace = TRUE, prob = proportion)
      RIi <- round(rlnorm(sampling, meanlog =  locations[paste0("mean", as.character(category))], sdlog = shapes[paste0("sd", as.character(category))])+1)
      rep_RIi <- c(rep_RIi, mean(1/RIi))
    }
  } else {
    rep_RIi <- NULL
    for (r in 1:replicates) {
      RIi <- sample(RI, sampling, replace = TRUE)
      rep_RIi <- c(rep_RIi, mean(1/RIi))
    }
  }
  return(c(mean=mean(rep_RIi), se=sd(rep_RIi)))
}


