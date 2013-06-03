#' plot.mcmcComposite plots the result of a MCMC search
#' @title Plot the result of a MCMC search
#' @author Marc Girondot
#' @return None
#' @param x A mcmcComposite object obtained after MHmcmc()
#' @param chain The chain to use
#' @param parameters Name of parameters or their number (see description)
#' @param ... Graphical parameters to be send to hist()
#' @description Plot the result of a MCMC search.\cr
#' The parameters to use can be called by:\cr
#' parameters="all"\cr
#' parameters=1:4\cr
#' parameters=c("PAR1", "PAR2", "PAR5")\cr
#' parameters=c(TRUE, TRUE, FALSE, TRUE)\cr
#' @examples 
#' library(phenology)
#' # Read a file with data
#' \dontrun{
#' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
#' }
#' data(Gratiot)
#' # Generate a formatted list named data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' \dontrun{
#' result_Gratiot<-fit_phenology(data=data_Gratiot, 
#' 		parametersfit=parg, parametersfixed=NULL, trace=1)
#' }
#' data(result_Gratiot)
#' # Generate set of priors for Bayesian analysis
#' \dontrun{
#' pmcmc <- phenology_MHmcmc_p(result_Gratiot)
#' }
#' pmcmc <- structure(c("dunif", "dunif", "dunif", "dunif", "dunif", "dunif", 
#' "dunif", "dunif", "0", "0", "0", "0", "0", "0", "0", "0", "200", 
#' "365", "200", "50", "200", "5", "5", "10", "2", "2", "2", "2", 
#' "2", "2", "2", "2", "0", "0", "0", "0", "0", "0", "0", "0", "200", 
#' "365", "200", "50", "200", "5", "5", "10", "95.826796339888", 
#' "175.36499338462", "62.4313052780003", "6.77668901451618e-05", 
#' "33.1138407661406", "0.21779065736816", "0.424368825094697", 
#' "3.58302217559733"), .Dim = c(8L, 7L), .Dimnames = list(c("LengthB", 
#' "Peak", "LengthE", "Flat", "Max_Gratiot", "MinB_Gratiot", "MinE_Gratiot", 
#' "Theta"), c("Density", "Prior1", "Prior2", "SDProp", "Min", "Max", 
#' "Init")))
#' \dontrun{
#' res_mcmc <- phenology_MHmcmc(result = result_Gratiot, n.iter = 10000, 
#' 		parametersMCMC = pmcmc, n.chains = 1, n.adapt = 0, thin = 1, trace = FALSE)
#' data(res_mcmc)
#' plot(res_mcmc, parameters=3, xlim=c(230, 330))
#' }
#' @method plot mcmcComposite
#' @export

plot.mcmcComposite<-function(x, ... , chain=1, parameters=1) {

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

	L <- modifyList(list(ylab="Density", xlab=rownames(Parameters)[[variable]], main="", freq=FALSE), list(x=x, ...)) 

	do.call(hist, L) 

	par(new=TRUE)

	x2 <- (par("usr")[1]+par("usr")[2]*26)/27
	x1 <- x2*26-par("usr")[2]/0.04

	plot(seq(from=x1, to=x2, length=100), get(Parameters[variable, "Density"])(seq(from=x1, to=x2, length=100), 
	as.numeric(Parameters[variable,"Prior1"]),as.numeric(Parameters[variable,"Prior2"])), type="l", col="red", xlab="", ylab="", axes = FALSE, xlim=c(x1, x2), main="")

	legend("topright", c("Prior", "Posterior"), lty=1, col=c('red', 'black'), bty = "n")

    }
  }
}
