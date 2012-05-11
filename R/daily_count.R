#' daily_count estimates nest number based on set of parameters.
#' @title Estimate expected counts based on set of parameters.
#' @author Marc Girondot
#' @param d Ordinal date
#' @param xpar Set of fixed+fitted parameters
#' @param print If TRUE, the result is displayed
#' @param help If TRUE, an help is displayed
#' @return The number of each day in d
#' @description Function estimates counts based on set of parameters.


.daily_count <-
function(d, xpar, print=TRUE, help=FALSE) {

.phenology.env<- NULL
rm(.phenology.env)


if(help) {
	cat("This function is for internal use only. The syntaxe is :\n")
	cat(".daily_count(day, xpar)\n")
	cat("With xpar<-c(par, pfixed) being the list of parameters.\n")
	cat("It must be preceeded by a call to xpar<-.format_par(xpar, serie, help=FALSE).\n")
} else {

	
# print(xpar)

nn<-ifelse(d<xpar["Begin"], xpar["MinB"],
		ifelse(d<xpar["PmoinsF"], ((1+cos(pi*(xpar["PmoinsF"]-d)/xpar["PmoinsFB"]))/2)*xpar["MaxMinB"]+xpar["MinB"],
			ifelse(d<xpar["PplusF"], xpar["Max"],
				ifelse(d<xpar["End"], ((1+cos(pi*(d-(xpar["PplusF"]))/xpar["EPplusF"]))/2)*xpar["MaxMinE"]+xpar["MinE"],
					xpar["MinE"]
				)
			)
		)
	)

if (any(is.na(nn))) {
	print("Global: Error, the parameters at the time of error are:")
	print(xpar)
	assign("par_error", xpar, envir=as.environment(.phenology.env))
	pause <- scan() 
}


if (xpar["sin"]) {
	ns<-sin(2*pi*((d+xpar["Delta"])/xpar["Phi"]))*(xpar["Alpha"]+(xpar["Beta"]*nn^xpar["Tau"]))
if (any(is.na(ns))) {
	print(d)
	print("Sin: Error, the parameters at the time of error are:")
	print(xpar)
	assign("par_error", xpar, envir=as.environment(.phenology.env))
	pause <- scan() 
}
} else {
	ns<-0
}

if (xpar["sin1"]) {
	ns1<-sin(2*pi*((d+xpar["Delta1"])/xpar["Phi1"]))*(xpar["Alpha1"]+(xpar["Beta1"]*nn^xpar["Tau1"]))
if (any(is.na(ns1))) {
	print(d)
	print("Sin 1: Error, the parameters at the time of error are:")
	print(xpar)
	assign("par_error", xpar, envir=as.environment(.phenology.env))
	pause <- scan() 
}
} else {
	ns1<-0
}


if (xpar["sin2"]) {
	ns2<-sin(2*pi*((d+xpar["Delta2"])/xpar["Phi2"]))*(xpar["Alpha2"]+(xpar["Beta2"]*nn^xpar["Tau2"]))
if (any(is.na(ns2))) {
	print(d)
	print("Sin 2: Error, the parameters at the time of error are:")
	print(xpar)
	assign("par_error", xpar, envir=as.environment(.phenology.env))
	pause <- scan() 
}
} else {
	ns2<-0
}

nn<-nn+ns+ns1+ns2
nn[is.na(nn)]<-1E-9

nn<-ifelse((nn<=0) & (d<xpar["Begin"]), xpar["MinB"],
		ifelse((nn<=0) & (d>xpar["Begin"]), xpar["MinE"],
			ifelse(nn<=0,(xpar["MinB"]+xpar["MinE"])/2,
				nn
			)
		)
	)

nn[nn==0]<-1E-9


# je suis en en mode interactif, j'affiche le résultat
if (print) {
	print(paste("Day ", d, "Number ", as.numeric(nn)))
}

	return(as.numeric(nn))
}
}

