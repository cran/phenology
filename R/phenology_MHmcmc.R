#' phenology_MHmcmc runs the Metropolis-Hastings algorithm for data (Bayesian MCMC)
#' @title Run the Metropolis-Hastings algorithm for data
#' @author Marc Girondot
#' @return A list with resultMCMC being mcmc.list object, resultLnL being likelihoods and parametersMCMC being the parameters used
#' @param n.iter Number of iterations for each step
#' @param parametersMCMC A set of parameters used as initial point for searching with information on priors
#' @param result An object obtained after a SearchR fit
#' @param n.chains Number of replicates
#' @param n.adapt Number of iterations before to store outputs
#' @param thin Number of iterations between each stored output
#' @param adaptive Should an adaptive process for SDProp be used
#' @param adaptive.lag  Lag to analyze the SDProp value in an adaptive content
#' @param adaptive.fun Function used to change the SDProp
#' @param trace TRUE or FALSE or period, shows progress
#' @param traceML TRUE or FALSE to show ML
#' @param filename If intermediate is not NULL, save intermediate result in this file
#' @param intermediate Period for saving intermediate result, NULL for no save
#' @param previous Previous result to be continued. Can be the filename in which intermediate results are saved.
#' @param WAIC Should WAIC information be saved?
#' @param WAIC.bybeach Should the WAIC be saved by beach or by count?
#' @description Run the Metropolis-Hastings algorithm for data.\cr
#' The number of iterations is n.iter+n.adapt+1 because the initial likelihood is also displayed.\cr
#' I recommend thin=10.\cr
#' If initial point is maximum likelihood, n.adapt = 0 is a good solution.\cr
#' The parameters intermediate and filename are used to save intermediate results every 'intermediate' iterations (for example 1000). Results are saved in a file of name filename.\cr
#' The parameter previous is used to indicate the list that has been save using the parameters intermediate and filename. It permits to continue a mcmc search.\cr
#' These options are used to prevent the consequences of computer crash or if the run is very very long and computer processes are time limited.
#' @family Phenology model
#' @examples 
#' \dontrun{
#' library(phenology)
#' data(Gratiot)
#' # Generate a formatted list named data_Gratiot 
#' data_Gratiot <- add_phenology(Gratiot, name="Complete", 
#'     reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg <- par_init(data_Gratiot, fixed.parameters=NULL)
#' # Run the optimisation
#' result_Gratiot <- fit_phenology(data=data_Gratiot, 
#' 		fitted.parameters=parg, fixed.parameters=NULL)
#' # Generate set of priors for Bayesian analysis
#' pmcmc <- phenology_MHmcmc_p(result_Gratiot, accept = TRUE)
#' # Here no WAIC data
#' result_Gratiot_mcmc <- phenology_MHmcmc(result = result_Gratiot, n.iter = 10000, 
#'                                         parametersMCMC = pmcmc, n.chains = 1, 
#'                                         n.adapt = 0, thin = 10, trace = FALSE, 
#'                                         WAIC=FALSE)
#' # Here the WAIC is at the level of the beaches
#' result_Gratiot_mcmc <- phenology_MHmcmc(result = result_Gratiot, n.iter = 10000, 
#'                                         parametersMCMC = pmcmc, n.chains = 1, 
#'                                         n.adapt = 0, thin = 10, trace = FALSE, 
#'                                         WAIC=TRUE, WAIC.bybeach=TRUE)
#' # Here the WAIC is at the level of the counts
#' result_Gratiot_mcmc <- phenology_MHmcmc(result = result_Gratiot, n.iter = 10000, 
#'                                         parametersMCMC = pmcmc, n.chains = 1, 
#'                                         n.adapt = 0, thin = 10, trace = FALSE, 
#'                                         WAIC=TRUE, WAIC.bybeach=FALSE)
#' loo::loo(result_Gratiot_mcmc$WAIC)
#' loo::waic(result_Gratiot_mcmc$WAIC)
#' # Get standard error of parameters
#' summary(result_Gratiot_mcmc)
#' # Make diagnostics of the mcmc results using coda package
#' mcmc <- as.mcmc(result_Gratiot_mcmc)
#' require(coda)
#' heidel.diag(mcmc)
#' raftery.diag(mcmc)
#' autocorr.diag(mcmc)
#' acf(mcmc[[1]][,"LengthB"], lag.max=200, bty="n", las=1)
#' acf(mcmc[[1]][,"Max_Complete"], lag.max=50, bty="n", las=1)
#' batchSE(mcmc, batchSize=100)
#' # The batch standard error procedure is usually thought to 
#' # be not as accurate as the time series methods used in summary
#' summary(mcmc)$statistics[,"Time-series SE"]
#' plot(result_Gratiot_mcmc, parameters="Max_Complete", las=1, xlim=c(-10, 210))
#' plot(result_Gratiot_mcmc, what="MarkovChain", ylim=c(20, 40))
#' plot(result_Gratiot_mcmc, what="LnL")
#' }
#' @export


phenology_MHmcmc<-function(result=stop("An output from fit_phenology() must be provided"), 
                           n.iter=10000, 
                           parametersMCMC=stop("A model generated with phenology_MHmcmc_p() must be provided"), 
                           n.chains = 1, 
                           n.adapt = 1000, 
                           thin=1, 
                           WAIC=FALSE, 
                           WAIC.bybeach=TRUE, 
                           trace=FALSE, 
                           traceML=FALSE, 
                           adaptive=TRUE, 
                           adaptive.lag=500, 
                           adaptive.fun=function(x) {ifelse(x>0.234, 1.3, 0.7)},
                           intermediate=NULL, 
                           filename="intermediate.Rdata", 
                           previous=NULL) {
  
  # result <- NULL; n.iter <- 10000; parametersMCMC <- NULL; n.chains = 1; n.adapt = 0; thin = 1; WAIC=FALSE; trace = FALSE; traceML = FALSE ; adaptive=FALSE; adaptive.lag=500; adaptive.fun=function(x) {ifelse(x>0.234, 1.3, 0.7)}; intermediate=NULL; filename="intermediate.Rdata"; previous=NULL
  # result <- result_Gratiot; parametersMCMC <- phenology_MHmcmc_p(result_Gratiot, accept = TRUE)
  
  mc.cores <- getOption("mc.cores", detectCores())
  
  message(paste0("I will use ", as.character(mc.cores)," cores."))
  
  
  if (is.character(previous)) {
    itr <- NULL
    load(previous)
    previous <- itr
    rm(itr)
    print("Continue previous mcmc run")
  }
  
  if (!inherits(result, "phenology")) {
    stop("An output of fit_phenology() must be provided")
  }
  
  pt <- list(data=result$data, fixed=result$fixed.parameters, 
             out=TRUE, 
             model_before=result$model_before, 
             cofactors=result$cofactors,
             add.cofactors=result$add.cofactors,
             zero=result$zero, 
             method_Snbinom=result$method_Snbinom, 
             store.intermediate=FALSE, 
             file.intermediate="", 
             WAIC=WAIC, 
             WAIC.bybeach=WAIC.bybeach)
  
  print(parametersMCMC)
  
  if (WAIC.bybeach) {
    n.datapoints <- length(result$data)
  } else {
    n.datapoints <- sum(unlist(lapply(result$data, FUN=function(x) nrow(x))))
  }
  
  out <- MHalgoGen(n.iter=n.iter, parameters=parametersMCMC, 
                   n.chains = n.chains, n.adapt = n.adapt, thin=thin, 
                   WAIC.out=WAIC, 
                   n.datapoints=n.datapoints, 
                   adaptive = adaptive, adaptive.fun = adaptive.fun, adaptive.lag = adaptive.lag,
                   trace=trace, traceML = traceML, pt=pt, 
                   likelihood=getFromNamespace(".Lnegbin", ns="phenology"))
  
  fin <- try(summary(out), silent=TRUE)
  
  if (inherits(fin, "try-error")) {
    lp <- rep(NA, nrow(out$parametersMCMC$parameters))
    names(lp) <- rownames(out$parametersMCMC$parameters)
    out <- c(out, SD=list(lp))
  } else {
    out <- c(out, SD=list(fin$statistics[,"SD"]))
  }
  
  out <-addS3Class(out, "mcmcComposite")
  
  return(out)
  
}
