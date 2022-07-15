#' ExponentialRegression is used to fit additive, multiplicative or mixte exponential regression
#' @title Non-biased exponential regression
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Return a list with the results of exponential regression
#' @param data A data.frame with a column time and a column number.
#' @param colname.time Name of the column to be used as time index.
#' @param colname.number Name of the column to be used as number index.
#' @param fitted.parameters A named vector with the parameters to be fitted
#' @param fixed.parameters  A named vector with the parameters to be fixed
#' @param weights	an optional vector of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. 
#' @param n.iter The number of MCMC iterations.
#' @param prior.default.distribution The default prior distribution; see description.
#' @param density by default is dnorm but can be dnbinom_new
#' @family Exponential regression
#' @description The idea of this function is to fit a regression exponential model using MCMC 
#' because regression model using glm can produce biased outputs.\cr
#' prior.default.distribution can be "dnorm", "dunif", or "dgamma". Note that if you propose 
#' dgamma prior, it will use uniform prior for r because r can be negative.\cr
#' The SD model is asd*Nt+bsd. The asd parameter represents multiplicative model and the 
#' bsd parameter represents additive model. Both can be used simultaneously.\cr
#' density can be dnorm or dnbinom_new (from HelpersMG package). dnbinom_new() is a 
#' negative binomial with mean and sd parametrization.
#' @examples
#' \dontrun{
#' library("phenology")
#' t <- 1:100
#' N0 <- 100
#' r <- 0.05
#' y <- N0*exp(t*r)
#' 
#' # Multiplicative model
#' Nt <- rnorm(100, mean=y, sd=0.2*y)
#' df <- data.frame(time=t, numbers=Nt)
#' g <- ExponentialRegression(data=df, fitted.parameters=c(N0=NA, r=NA, asd=NA))
#' plot(g, parameters="r")
#' as.parameters(g, index="median")
#' # Note that if you propose gamma prior, it will use uniform prior for r
#' # because r can be negative
#' g <- ExponentialRegression(data=df, 
#'                            fitted.parameters=c(N0=NA, r=NA, asd=NA), 
#'                            prior.default.distribution="dgamma")
#' plot(g, parameters="r")
#' as.parameters(g, index="median")
#' 
#' # Additive model
#' Nt <- rnorm(100, mean=y, sd=5)
#' df <- data.frame(time=t, numbers=Nt)
#' g <- ExponentialRegression(data=df, fitted.parameters=c(N0=NA, r=NA, bsd=NA))
#' plot(g, parameters="r")
#' as.parameters(g, index="median")
#' 
#' # Mixt model
#' Nt <- rnorm(100, mean=y, sd=0.2*y+5)
#' df <- data.frame(time=t, numbers=Nt)
#' g <- ExponentialRegression(data=df, fitted.parameters=c(N0=NA, r=NA, asd=NA, bsd=NA))
#' plot(g, parameters="r")
#' as.parameters(g, index="median")
#' 
#' # Example with 3 common ways to perform the regression
#' t <- 1:100
#' N0 <- 100
#' r <- 0.05
#' y <- N0*exp(t*r)
#' out_glm <- NULL
#' out_mcmc <- NULL
#' out_nls <- NULL
#' for (i in 1:500) {
#'         print(i)
#'         set.seed(i)
#'         Nt <- rnorm(100, mean=y, sd=0.2*y)
#'         df <- data.frame(time=t, numbers=Nt)
#'         g0 <- glm(log(numbers) ~ time, data = df)
#'         out_glm <- c(out_glm, c(exp(coef(g0)[1]), coef(g0)[2]))
#'         g1 <- ExponentialRegression(data=df, n.iter=c(10000, 20000))
#'         out_mcmc <- c(out_mcmc, as.parameters(g1, index="median")[1:2])
#'         g2 <- nls(numbers ~ N0*exp(r*time), start = list(N0 = 100, r = 0.05), data = df)
#'         out_nls <- c(out_nls, coef(g2))
#' }
#' # In conclusion the method proposed here has no biais as compare to glm and nls fits
#' out_glm <- matrix(out_glm, ncol=2, byrow=TRUE)
#' out_mcmc <- matrix(out_mcmc, ncol=2, byrow=TRUE)
#' out_nls <- matrix(out_nls, ncol=2, byrow=TRUE)
#' mean(out_glm[, 1]); mean(out_mcmc[, 1]); mean(out_nls[, 1])
#' sd(out_glm[, 1])/sqrt(nrow(out_glm)); sd(out_mcmc[, 1])/sqrt(nrow(out_mcmc)); 
#' sd(out_nls[, 1])/sqrt(nrow(out_nls))
#' mean(out_glm[, 2]); mean(out_mcmc[, 2]); mean(out_nls[, 2])
#' sd(out_glm[, 2])/sqrt(nrow(out_glm)); sd(out_mcmc[, 2])/sqrt(nrow(out_mcmc)); 
#' sd(out_nls[, 2])/sqrt(nrow(out_nls))
#' }
#' @export

ExponentialRegression <- function (data=stop("A data.frame with values"), 
                                   colname.time="time", 
                                   colname.number="numbers", 
                                   weights=NULL, 
                                   fitted.parameters=c(N0=NA, r=NA, asd=NA), 
                                   fixed.parameters=NULL, 
                                   n.iter=c(100000, 100000), 
                                   prior.default.distribution="dnorm", 
                                   density=dnorm) {
  
  # colname.time="time"; colname.number="numbers"; weights=NULL; fitted.parameters=c(N0=NA, r=NA, asd=NA); fixed.parameters=NULL; n.iter=100000; prior.default.distribution="dnorm"
  
  if (length(n.iter) == 1) n.iter <- rep(n.iter, 2)
  
  if (is.null(weights)) weights <- rep(1, nrow(data)) 
  if (length(weights) != nrow(data)) {
    warning("The weights are forced to be the same for all data.")
    weights <- rep(1, nrow(data)) 
  }
  colnames(data)[which(colnames(data)==colname.time)] <- "time"
  colnames(data)[which(colnames(data)==colname.number)] <- "numbers"
  
  lnexp <- getFromNamespace(".lnexp", ns="phenology")
  
  g0 <- glm(log(numbers) ~ time, data = data, weights = weights)
  
  low <- NULL
  high <- NULL
  fit <- NULL
  rules <- NULL
  
  if (any(names(fitted.parameters) == "N0")) {
    fit <- c(fit, N0=unname(exp(coef(g0)["(Intercept)"])))
    low <- c(low, N0=1E-9)
    high <- c(high, N0=unname(exp(coef(g0)["(Intercept)"]*2)))
    rules <- rbind(rules, 
                   data.frame(Name="N0", Min=1E-9, 
                              Max=unname(exp(coef(g0)["(Intercept)"]*2))))
  }
  
  if (any(names(fitted.parameters) == "r")) {
    fit <- c(fit, r=unname(abs(coef(g0)["time"])))
    low <- c(low, r=unname(-abs(coef(g0)["time"]*2)))
    high <- c(high, r=unname(abs(coef(g0)["time"]*2)))
    rules <- rbind(rules, 
                   data.frame(Name="r", Min=unname(-abs(coef(g0)["time"]*2)), 
                              Max=unname(abs(coef(g0)["time"]*2))))
  }
  
  if (any(names(fitted.parameters) == "asd")) {
    fit <- c(fit, asd=summary(g0)$dispersion)
    low <- c(low, asd=1E-9)
    high <- c(high, asd=summary(g0)$dispersion*2)
    rules <- rbind(rules, 
                   data.frame(Name="asd", Min=1E-9, 
                              Max=summary(g0)$dispersion*2))
    
  }
  
  if (any(names(fitted.parameters) == "bsd")) {
    fit <- c(fit, bsd=unname(exp(coef(g0)["(Intercept)"]))/2)
    low <- c(low, bsd=1E-9)
    high <- c(high, bsd=unname(exp(coef(g0)["(Intercept)"])))
    rules <- rbind(rules, 
                   data.frame(Name="bsd", Min=1E-9, 
                              Max=unname(exp(coef(g0)["(Intercept)"]))))
    
  }
  
  # lnexp(df=data, par=fit, fixed.parameters=fixed.parameters)
  
  G <- optim(par=fit, 
             fixed.parameters=fixed.parameters, 
             method = "L-BFGS-B", 
             density = density, 
             lower = low, 
             upper = high,
             df=data, weights=weights, fn=lnexp, hessian = TRUE)
  
  prior <- setPriors(par = G$par, 
                     se = SEfromHessian(a=G$hessian), 
                     density = prior.default.distribution, rules = rules, silent = TRUE)
  
  prior$SDProp <- G$par/10
  
  MCMC <- MHalgoGen(likelihood = lnexp, parameters_name = "par", 
                    fixed.parameters=fixed.parameters, 
                    density = density, 
                    parameters = prior, df=data, adaptive = TRUE, n.iter = n.iter[1])
  prior$SDProp <- MCMC$parametersMCMC$SDProp.end
  prior$Init <- as.parameters(MCMC, index = "median")
  MCMC <- MHalgoGen(likelihood = lnexp, parameters_name = "par", 
                    fixed.parameters=fixed.parameters, 
                    density = density, 
                    parameters = prior, df=data, adaptive = FALSE, n.iter = n.iter[2])
  return(MCMC)
}

.lnexp <- function(par, fixed.parameters=NULL, df, weights=1, density=dnorm) {
  par <- c(par, fixed.parameters)
  if (is.na(par["asd"])) par <- c(par, asd=0)
  if (is.na(par["bsd"])) par <- c(par, bsd=0)
  y <- par["N0"]*exp(df[, "time"]*par["r"])
  # Not clear if I should use sd=df[, "y"] or sd=par["sd"]*y
  # I use y because it is always positif
  LnL <- -sum(density(df[, "numbers"], mean=y, sd=par["asd"]*y+par["bsd"], log = TRUE)*weights)
  return(LnL)
}
