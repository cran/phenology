#' phenology_MHmcmc_p generates set of parameters to be used with MHmcmc()
#' @title Generates set of parameters to be used with phenology_MHmcmc()
#' @author Marc Girondot
#' @return A matrix with the parameters
#' @param result An object obtained after a SearchR fit
#' @description Interactive script used to generate set of parameters to be used with phenology_MHmcmc().\cr
#' @examples 
#' library(phenology)
#' # Read a file with data
#' # Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", , header=FALSE)
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot<-add_format(origin=NULL, add=Gratiot, name="Complete", reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' ## not run
#' # result_Gratiot<-fit_phenology(data=data_Gratiot, parametersfit=parg, parametersfixed=NULL, trace=1)
#' ## end not run
#' data(result_Gratiot)
#' # Generate set of priors for Bayesian analysis
#' ## not run
#' ## pmcmc <- phenology_MHmcmc_p(result_Gratiot)
#' ## end not run
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
#' @export

phenology_MHmcmc_p<-function(result=stop("An output from fit_phenology() must be provided")) {

if (class(result)!="phenology") {
	cat("An output from fit_phenology() must be provided\n")
	return()
}

# d'abord je sors les paramètres à utiliser

par <- result$par

# "Peak"
Peak <- c("dunif", 0, 365, 2, 0, 365, ifelse(is.na(par["Peak"]), 160, par["Peak"]))

# "Flat"
Flat <- c("dunif", 0, 50, 2, 0, 50, ifelse(is.na(par["Flat"]), 5, par["Flat"]))

# "Begin"
Begin <- c("dunif", 0, 365, 2, 0, 365, ifelse(is.na(par["Begin"]), 160, par["Begin"]))

# "End"
End <- c("dunif", 0, 365, 2, 0, 365, ifelse(is.na(par["End"]), 220, par["End"]))

# "Length"
Length <- c("dunif", 0, 200, 2, 0, 200, ifelse(is.na(par["Length"]), 100, par["Length"]))

# "LengthE"
LengthE <- c("dunif", 0, 200, 2, 0, 200, ifelse(is.na(par["LengthE"]), 100, par["LengthE"]))

# "LengthB"
LengthB <- c("dunif", 0, 200, 2, 0, 200, ifelse(is.na(par["LengthB"]), 100, par["LengthB"]))

# "Max"
Max <- c("dunif", 0, 200, 2, 0, 200, ifelse(is.na(par["Max"]), 100, par["Max"]))

# "PMin"
PMin <- c("dunif", 0, 10, 2, 0, 10, ifelse(is.na(par["PMin"]), 2, par["Pmin"]))

# "Min"
Min <- c("dunif", 0, 5, 2, 0, 5, ifelse(is.na(par["Min"]), 1, par["Min"]))

# "PMinE"
PMinE <- c("dunif", 0, 10, 2, 0, 10, ifelse(is.na(par["PMinE"]), 2, par["PminE"]))

# "MinE"
MinE <- c("dunif", 0, 5, 2, 0, 5, ifelse(is.na(par["MinE"]), 1, par["MinE"]))

# "PMinB"
PMinB <- c("dunif", 0, 10, 2, 0, 10, ifelse(is.na(par["PMinB"]), 2, par["PminB"]))

# "MinB"
MinB <- c("dunif", 0, 5, 2, 0, 5, ifelse(is.na(par["MinB"]), 1, par["MinB"]))

# "Phi"
Phi <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Phi"]), 10, par["Phi"]))

# "Delta"
Delta <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Delta"]), 10, par["Delta"]))

# "Alpha"
Alpha <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Alpha"]), 10, par["Alpha"]))

# "Beta"
Beta <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Beta"]), 10, par["Beta"]))

# "Tau"
Tau <- c("dunif", 0, 2, 2, 0, 2, ifelse(is.na(par["Tau"]), 1, par["Tau"]))

# "Phi1"
Phi1 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Phi1"]), 10, par["Phi1"]))

# "Delta1"
Delta1 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Delta1"]), 10, par["Delta1"]))

# "Alpha1"
Alpha1 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Alpha1"]), 10, par["Alpha1"]))

# "Beta1"
Beta1 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Beta1"]), 10, par["Beta1"]))

# "Tau1"
Tau1 <- c("dunif", 0, 2, 2, 0, 2, ifelse(is.na(par["Tau1"]), 1, par["Tau1"]))

# "Phi2"
Phi2 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Phi2"]), 10, par["Phi2"]))

# "Delta2"
Delta2 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Delta2"]), 10, par["Delta2"]))

# "Alpha2"
Alpha2 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Alpha2"]), 10, par["Alpha2"]))

# "Beta2"
Beta2 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Beta2"]), 10, par["Beta2"]))

# "Tau2"
Tau2 <- c("dunif", 0, 2, 2, 0, 2, ifelse(is.na(par["Tau2"]), 1, par["Tau2"]))

# "Theta"
Theta <- c("dunif", 0, 10, 2, 0, 10, ifelse(is.na(par["Theta"]), 10, par["Theta"]))


priors <- list(Peak, Flat, Begin, End, Length, LengthE, LengthB, 
Length, PMin, Min, PMinE, MinE, PMinB, MinB, Phi, Delta, Alpha, 
Beta, Tau, Phi1, Delta1, Alpha1, Beta1, Tau1, Phi2, Delta2, 
Alpha2, Beta2, Tau2, Theta)

names(priors) <- c("Peak", "Flat", "Begin", "End", "Length", 
"LengthE", "LengthB", "Length", "PMin", "Min", "PMinE", "MinE", 
"PMinB", "MinB", "Phi", "Delta", "Alpha", "Beta", "Tau", "Phi1", 
"Delta1", "Alpha1", "Beta1", "Tau1", "Phi2", "Delta2", "Alpha2", 
"Beta2", "Tau2", "Theta")

for (i in 1:length(par)) {

		if (substr(names(par[i]), 1, 4)=="Max_") {
			priors <- c(priors, list(c("dunif", 0, 200, 2, 0, 200, par[i])))
			names(priors)[length(priors)] <- names(par[i])
		}
		
		
		if (substr(names(par[i]), 1, 4)=="Min_") {
			priors <- c(priors, list(c("dunif", 0, 5, 2, 0, 5, par[i])))
			names(priors)[length(priors)] <- names(par[i])
		}
		
		if (substr(names(par[i]), 1, 5)=="MinE_") {
			priors <- c(priors, list(c("dunif", 0, 5, 2, 0, 5, par[i])))
			names(priors)[length(priors)] <- names(par[i])
		}
		
		if (substr(names(par[i]), 1, 5)=="MinB_") {
			priors <- c(priors, list(c("dunif", 0, 5, 2, 0, 5, par[i])))
			names(priors)[length(priors)] <- names(par[i])
		}
		if (substr(names(par[i]), 1, 5)=="Peak_") {
			priors <- c(priors, list(c("dunif", 0, 365, 2, 0, 365, par[i])))
			names(priors)[length(priors)] <- names(par[i])
		}
}


prencours <- NULL

for (i in 1:length(par)) {

	prencours <- c(prencours, priors[[names(par)[i]]])
}



parametersMCMC <- matrix(prencours, ncol=7, byrow=T)
colnames(parametersMCMC) <- c("Density", "Prior1", "Prior2", "SDProp", "Min", "Max", "Init")
rownames(parametersMCMC)<-names(par)



parameters <- parametersMCMC

	repeat {

cat("Proposition:\n")
print(parameters)
cat("Name of the parameter to change or Enter to quit:\n")
f<-scan(nmax=1, quiet=TRUE, what=character())

if (length(f)==0) f <- "q"

if (f=="q") {
	return(parameters)
	
} else {

	variable <- which(f==names(par))
	if (length(variable)==0) {
	cat("The parameter does not exist:\n")
	} else {
	print(variable)
	cat(paste("Change for the parameter ",names(par)[variable],":\n",sep=""))

	cat(paste("Distribution of the prior (Enter for default ",parameters[variable, "Density"], "):", sep=""))
	density<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(density)!=0) { parameters[variable, "Density"] <- density } else { density <- parameters[variable, "Density"] }
	
	if (density == "dunif") {
	
	cat(paste("Distribution of the prior, Minimum (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior1"] <- f
	cat(paste("Distribution of the prior, Maximum (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior2"] <- f
	
	} else {
	
	if (density == "dnorm") {
	
	cat(paste("Distribution of the prior, Mean (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior1"] <- f
	cat(paste("Distribution of the prior, Standard deviation (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior2"] <- f
	
	} else {

	cat(paste("Distribution of the prior, value 1 (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior1"] <- f
	cat(paste("Distribution of the prior, value 2 (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior2"] <- f

	}
	}
	
	
	cat(paste("SD of new proposition (Enter for default ",parameters[variable, "SDProp"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "SDProp"] <- f
	cat(paste("Minimum for the parameter (default ",parameters[variable, "Min"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Min"] <- f
	cat(paste("Maximum for the parameter (Enter for default ",parameters[variable, "Max"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Max"] <- f
	cat(paste("Initial value (Enter for default ",parameters[variable, "Init"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Init"] <- f
	}

}


}

}