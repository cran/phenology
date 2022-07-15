#' Bayesian.remigration fits a remigration interval using Bayesian MCMC
#' @title Return a posterior remigration interval.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Return a posterior remigration interval.\cr
#' @param parameters Priors for Bayesian MCMC
#' @param data Data to be fitted
#' @param kl Maximum number of years for remigration intervals.
#' @param n.iter Number of iterations for MCMC
#' @param n.chains Number of replicates
#' @param n.adapt Number of iterations before to store outputs
#' @param thin Number of iterations between each stored output
#' @param adaptive Should an adaptive process for SDProp be used
#' @param adaptive.lag  Lag to analyze the SDProp value in an adaptive context
#' @param adaptive.fun Function used to change the SDProp
#' @param trace Or FALSE or period to show progress
#' @param filename If intermediate is not NULL, save intermediate result in this file
#' @param intermediate Period for saving intermediate result, NULL for no save
#' @param previous Previous result to be continued. Can be the filename in which intermediate results are saved.
#' @description Model of remigration interval\cr
#' @family Model of Remigration Interval
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' 
#' # Each year a fraction of 0.9 is surviving
#' s <- c(s=0.9)
#' # Probability of tag retention; 0.8
#' t <- c(t=0.8)
#' # Time-conditional return probability - This is the true remigration rate
#' r <- c(r1=0.1, r2=0.8, r3=0.7, r4=0.7, r5=1)
#' # Capture probability
#' p <- c(p1=0.6, p2=0.6, p3=0.6, p4=0.6, p5=0.6)
#' # Number of observations for 400 tagged females after 1, 2, 3, 4, and 5 years
#' OBS <- c(400, 10, 120, 40, 20, 10)
#' 
#' kl_s <- length(s)
#' kl_t <- length(t)
#' kl_r <- length(r)
#' kl_p <- length(p)
#' 
#' pMCMC <- data.frame(Density=c("newdbeta", "newdbeta", rep("dunif", kl_r), 
#'                               rep("newdbeta", kl_p), "dunif"), 
#'                     Prior1=c(s, t, rep(0, kl_r), rep(0.2, kl_p), 0), 
#'                     Prior2=c(0.02, 0.02, rep(1, kl_r), rep(0.08, kl_p), 10), 
#'                     SDProp=c(0.05, 0.05, rep(0.05, kl_r), rep(0.05, kl_p), 0.05), 
#'                  Min=c(0, 0, rep(0, kl_r), rep(0, kl_p), 0),  
#'                  Max=c(1, 1, rep(1, kl_r), rep(1, kl_p), 10),  
#'                  Init=c(s, t, r, p, 1), stringsAsFactors = FALSE, 
#'                  row.names=c("s", 
#'                                 "t", 
#'                                 names(r), 
#'                                 names(p), "sd")
#' )
#' rMCMC <- Bayesian.remigration(parameters = pMCMC, 
#' n.iter = 1000000, 
#' n.adapt = 300000,
#' trace=10000, 
#' data=OBS)
#' 
#' plot(rMCMC)
#' 
#' }
#' @export



Bayesian.remigration <- function(parameters = stop("Priors must be supplied"), 
                                 data=stop("data must be supplied"), 
                                 kl=NULL, 
                                 n.iter = 100000,
          n.chains = 1, n.adapt = 10000, thin = 1, trace = 10,
          adaptive = TRUE, adaptive.lag = 500, adaptive.fun = function(x) {    
            ifelse(x > 0.234, 1.3, 0.7) }, intermediate = NULL,
          filename = "intermediate.Rdata", previous = NULL) {
  
  rMCMC <- MHalgoGen(likelihood = LnRI_norm, 
                     data=data, 
                     kl=kl, 
                     parameters = parameters, 
                     n.adapt = n.adapt, 
                     n.iter = n.iter, 
                     trace=trace, 
                     adaptive = adaptive, 
                     thin = thin, 
                     n.chains = n.chains, 
                     adaptive.lag = adaptive.lag, 
                     adaptive.fun = adaptive.fun, 
                     intermediate = intermediate,
                     filename = filename, 
                     previous = previous)
  
  
  
  ri <- apply(rMCMC$resultMCMC[[1]][, substr(colnames(rMCMC$resultMCMC[[1]]), 1, 1)=="r"], 
              MARGIN = 1, FUN=function(r) {
                r <- r[order(as.numeric(substr(names(r), 2, nchar(r))))]
                rp <- rep(NA, length(r))
                rp[1] <- r[1]
                if (length(r) != 1) {
                  for (i in 2:length(r)) rp[i] <- r[i]*prod(1 - r[1:(i-1)])
                }
                sum(rp*(1:length(rp)))/sum(rp)
              })
  
  rMCMC <- list(rMCMC, RI=ri)
  rMCMC <- addS3Class(rMCMC, "Remigration")
  return(rMCMC)
}


