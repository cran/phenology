#' plot.mcmcComposite plots the result of a MCMC search
#' @title Plot the result of a MCMC search
#' @author Marc Girondot
#' @return None
#' @param x A mcmcComposite object obtained after phenology_fonctionMCMC()
#' @param chain The chain to use
#' @param parameters Name of parameters or their number (see description)
#' @param legend If FALSE, the legend is not shown
#' @param scale.prior If TRUE, the prior is scaled at the same size as posterior
#' @param ... Graphical parameters to be send to hist()
#' @description Plot the result of a MCMC search.\cr
#' The parameters to use can be called by:\cr
#' parameters="all"\cr
#' parameters=1:4\cr
#' parameters=c("PAR1", "PAR2", "PAR5")\cr
#' parameters=c(TRUE, TRUE, FALSE, TRUE)\cr
#' @examples 
#' \dontrun{
#' library(phenology)
#' data(Gratiot)
#' # Generate a formatted list named data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#'     reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' result_Gratiot<-fit_phenology(data=data_Gratiot, 
#' 		parametersfit=parg, parametersfixed=NULL, trace=1)
#' # Generate set of priors for Bayesian analysis
#' pmcmc <- phenology_MHmcmc_p(result_Gratiot, accept = TRUE)
#' result_Gratiot_mcmc <- phenology_MHmcmc(result = result_Gratiot, n.iter = 10000, 
#' parametersMCMC = pmcmc, n.chains = 1, n.adapt = 0, thin = 1, trace = FALSE)
#' # Get standard error of parameters
#' summary(result_Gratiot_mcmc)
#' # Make diagnostics of the mcmc results using coda package
#' mcmc <- as.mcmc(result_Gratiot_mcmc)
#' require(coda)
#' heidel.diag(mcmc)
#' raftery.diag(mcmc)
#' autocorr.diag(mcmc)
#' acf(mcmc[[1]][,"LengthB"], lag.max=200, bty="n", las=1)
#' acf(mcmc[[1]][,"Max_Gratiot"], lag.max=50, bty="n", las=1)
#' batchSE(mcmc, batchSize=100)
#' # The batch standard error procedure is usually thought to 
#' # be not as accurate as the time series methods used in summary
#' summary(mcmc)$statistics[,"Time-series SE"]
#' plot(result_Gratiot_mcmc, parameters=3, las=1, xlim=c(-10, 300))
#' }
#' @method plot mcmcComposite
#' @import coda
#' @export

plot.mcmcComposite <- function(x, ... , chain=1, parameters=1, scale.prior=FALSE, legend=TRUE) {

resultMCMC <- x

mcmc <- resultMCMC[["resultMCMC"]]
Parameters_L <- resultMCMC[["parametersMCMC"]]
Parameters <- Parameters_L[[1]]
n.iter <- Parameters_L[[2]]
n.chains <- Parameters_L[[3]]
n.adapt <- Parameters_L[[4]]
thin <- Parameters_L[[5]]

possible <- rownames(Parameters)
NbTS <- length(possible)

if (parameters[[1]]=="all") {
	series<-rep(TRUE, NbTS)
} else {
	if (any(!is.logical(parameters))) {
		if (is.numeric(parameters)) {
# Même si un nombre de paramètres plus important est indiqué, ne provoque pas d'erreur
			series <- rep(FALSE, max(NbTS, length(parameters)))
			series[parameters] <- TRUE
			series <- series[1:NbTS]
		} else {
			series <- (possible==parameters)
		}
	} else {
# c'est des valeurs logiques, je verifie si le bon nombre, sinon je recycle
		if (length(parameters)!=NbTS) {
			series <- rep(parameters, NbTS)
			series <- series[1:NbTS]
		}
	}
}

seriesx <- series

nbseries <- length(seriesx[seriesx==TRUE])

for(variable in 1:NbTS) {

if (seriesx[variable]) {

nitercorrige <- floor(n.iter/thin)
vals <- mcmc[[chain]][,variable]
# c'est quoi ça ??? 6/10/2012
# vals <- vals[(length(vals)-nitercorrige):length(vals)]

	x <- vals

if (Parameters[variable, "Density"]=="dunif") {
  # 22/8/2014
  # je suis en uniforme. Il faut que je limite les breaks
  
  mx <- as.numeric(Parameters[variable,"Prior2"])
  mn <- as.numeric(Parameters[variable,"Prior1"])
  mxlog <- 10^(floor(log10(mx))-1)
  xl <- as.numeric(c(mn, mx))+c(-mxlog, mxlog)
  br1 <- mn
  br2 <- mx
  interval <- (mx-mn)/20
  decalage <- br1 %% interval
  br <- seq(from=decalage, to=mx, by=interval)
  
	L <- modifyList(list(ylab="Density", xlab=rownames(Parameters)[[variable]], main="", freq=FALSE, 
	xlim=xl, breaks=br), list(x=x, ...)) 

  
} else {
	L <- modifyList(list(ylab="Density", xlab=rownames(Parameters)[[variable]], main="", freq=FALSE), list(x=x, ...)) 
}

	do.call(hist, L) 
	
	scl <- ScalePreviousPlot()

  sequence <- seq(from=scl$xlim[1], to=scl$xlim[2], length=200)
  p1 <- as.numeric(Parameters[variable,"Prior1"])
  p2 <- as.numeric(Parameters[variable,"Prior2"])
  y <- get(Parameters[variable, "Density"])(sequence, p1, p2)
  yl <- c(0, max(y))
  
if (Parameters[variable, "Density"]!="dunif") {
  par(new=TRUE)
  if (scale.prior) {
    plot(sequence, y, type="l", col="red", axes=FALSE, bty="n", xlab="", ylab="", main="", ylim=yl)
  } else {
    plot(sequence, y, type="l", col="red", axes=FALSE, bty="n", xlab="", ylab="", main="", ylim=scl$ylim[1:2])
  }
} else {
  if (scale.prior) {
  segments(x0=scl$xlim[1], x1=p1, y0=0, y1=0, col="red")
  segments(x0=p1, x1=p1, y0=0, y1=scl$ylim[2], col="red")
  segments(x0=p1, x1=p2, y0=scl$ylim[2], y1=scl$ylim[2], col="red")
  segments(x0=p2, x1=p2, y0=scl$ylim[2], y1=0, col="red")
  segments(x0=p2, x1=scl$xlim[2], y0=0, y1=0, col="red")
  } else {
    segments(x0=scl$xlim[1], x1=p1, y0=0, y1=0, col="red")
    segments(x0=p1, x1=p1, y0=0, y1=yl[2], col="red")
    segments(x0=p1, x1=p2, y0=yl[2], y1=yl[2], col="red")
    segments(x0=p2, x1=p2, y0=yl[2], y1=0, col="red")
    segments(x0=p2, x1=scl$xlim[2], y0=0, y1=0, col="red")
  }
}

if (legend) {
	legend("topright", c("Prior", "Posterior"), lty=1, col=c('red', 'black'), bty = "n")
}
    }
  }
}
