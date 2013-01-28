#' map_phenology generates a likelihood map.
#' @title Generate a likelihood map varying Phi and Delta.
#' @author Marc Girondot
#' @return Display a likelihood map
#' @param data dataset generated with add_format
#' @param parametersfixed Set of fixed parameters
#' @param parametersfit Set of parameters to be fitted
#' @param Phi Phi values to be analyzed
#' @param Delta Delta value to be analyzed
#' @param method_incertitude 2 [default] is the correct one from a statistical point of view; \cr
#'                           0 is an aproximate method more rapid; \cr
#'                           1 is an alternative more rapid but biased.\cr
#' @param zero_counts Example c(TRUE, TRUE, FALSE) indicates whether the zeros have 
#'                    been recorder for each of these timeseries. Defaut is TRUE for all.
#' @param progressbar If FALSE, do not show the progress bar
#' @param help If TRUE, an help is displayed
#' @description This function generates a map of likelihood varying Phi and Delta.\cr
#' When Delta is not given, the same precision as Phi is used.
#' @examples
#' library("phenology")
#' # Read a file with data
#' # Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", , header=FALSE)
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' # result_Gratiot<-fit_phenology(data=data_Gratiot, parametersfit=parg, parametersfixed=NULL, trace=1)
#' data(result_Gratiot)
#' # Extract the fitted parameters
#' parg1<-extract_result(result_Gratiot)
#' # Add constant Alpha and Tau values [day d amplitude=(Alpha+Nd*Beta)^Tau with Nd being the number of counts for day d]
#' pfixed<-c(parg1, Alpha=0, Tau=1)
#' pfixed<-pfixed[-which(names(pfixed)=="Theta")]
#' # The only fitted parameter will be Beta
#' parg2<-c(Beta=0.5, parg1["Theta"])
#' # Generate a likelihood map [default Phi=seq(from=0.1, to=20, length.out=100) but it is very long]
#' # Take care, it takes 20 hours ! The data map_Gratiot has the result
#' # map_Gratiot<-map_phenology(data=data_Gratiot, Phi=seq(from=0.1, to=20, length.out=100), parametersfit=parg2, parametersfixed=pfixed)
#' data(map_Gratiot)
#' # Plot the map
#' plot(map_Gratiot, col=heat.colors(128))
#' # Plot the min(-Ln L) for Phi varying at any delta value
#' plot_phi(map=map_Gratiot)
#' # Plot the min(-Ln L) for Delta varying with Phi equal to the value for maximum likelihood
#' plot_delta(map=map_Gratiot)
#' # Plot the min(-Ln L) for Delta varying with Phi the nearest to 15
#' plot_delta(map=map_Gratiot, Phi=15)
#' @export

map_phenology <-
function(data=NULL, parametersfit=NULL, parametersfixed=NA, Phi=seq(from=0.2,to=20, length.out=100), Delta=NULL, method_incertitude=2, zero_counts=TRUE, progressbar=TRUE, help=FALSE) {
if (help || is.null(data)) {
	cat("This function generates a map of likelihood varying Phi and Delta.\n")
	cat("map<-map_phenology(data=dataset, parametersfit=par, parametersfixed=pfixed,\n")
	cat("+      Phi=seq(from=0.2,to=20, length.out=100), Delta=seq(from=0, to=5, length.out=101),\n")	
	cat("+      method_incertitude=2, zero_counts=TRUE)\n")
	cat("or if no parameter is fixed and default is used:\n")
	cat("map<-map_phenology(data=dataset, parametersfit=par)\n")
	cat("or\n")	
	cat("map_phenology(help=TRUE) to have this help !\n")
	cat("parameters are the same than for the fit_phenology() function except for trace that is disabled.\n")
	cat("If Alpha, Beta or Tau are not indicated, Alpha and Tau are set to 0 and 1 and Beta is fitted.\n")
	cat("Only one set of Alpha, Beta, Tau, Phi and Delta are used for all timeseries present in data.\n")
	cat("Note that it is possible to fit or fixed Alpha[n], Beta[n], Tau[n], Phi[n] and Delta[n] with [n]=1 or 2\n")
	cat("and then it is possible to use this function to establish the likelihood map for a\n")
	cat("second or third sinusoids added to the global pattern.\n")
	cat("If Delta is not specified, it is estimated from Phi.\n")
	cat("See Girondot, M., Rivalan, P., Wongsopawiro, R., Briane, J.-P., Hulin, V., Caut, S., Guirlet, E. \n")
	cat("& Godfrey, M.H. (2006) Phenology of marine turtle nesting revealed by a statistical model of the \n")
	cat("nesting season. BMC Ecology, 6, 11.\n")
	cat("for an exemple.\n")

} else {

#.phenology.env<- NULL
#rm(.phenology.env)

if (is.null(parametersfixed)) {parametersfixed<-NA}
if (is.null(parametersfit)) {parametersfit<-NA}

# je varie Phi et Delta

#create 2 vectors in form of numeric sequence, for Delta and Phi
Phivalue=Phi
if (is.null(Delta)) {

Deltavalue=seq(from=0, to=max(Phivalue)/2, length.out=length(Phivalue)+1)

} else {
	Deltavalue=Delta
}

LPhi<-length(Phivalue)
LDelta<-length(Deltavalue)

# SET MATRIX
matrix(data=NA, LPhi, LDelta) ->input

#	.phenology.env <<- new.env(parent=.GlobalEnv)
#	assign("data", data, envir=.phenology.env)
#	assign("fixed", parametersfixed, envir=.phenology.env)
#	assign("incertitude", method_incertitude, envir=.phenology.env)


	if (length(zero_counts)==1) {zero_counts<-rep(zero_counts, length(data))}
	if (length(zero_counts)!=length(data)) {
		print("zero_counts parameter must be TRUE (the zeros are used for all timeseries) or FALSE (the zeros are not used for all timeseries) or possess the same number of logical values than the number of series analyzed.")
		return()
	}

#	assign("zerocounts", zero_counts, envir=.phenology.env)


# si ni Alpha ni Beta ne sont à ajuster, je mets Beta
if (is.na(parametersfit["Alpha"]) && is.na(parametersfit["Beta"])) {
	if (all(is.na(parametersfit))) {
		parametersfit<-c(Beta=0)
	} else {
		parametersfit<-c(parametersfit, Beta=0)
	}
}

# si Beta est à la fois fixe et à ajuster, je le retire des fixes
if (!is.na(parametersfit["Beta"]) && !is.na(parametersfixed["Beta"])) {parametersfixed<-parametersfixed[!names(parametersfixed)=="Beta"]}

# je vérifie que Alpha, Beta et Tau apparaissent bien au moins une fois, sinon je les mets en fixe
xpar<-c(parametersfit, parametersfixed)
if (is.na(xpar["Alpha"])) {parametersfixed<-c(parametersfixed, Alpha=0)}
if (is.na(xpar["Beta"])) {parametersfixed<-c(parametersfixed, Beta=0)}
if (is.na(xpar["Tau"])) {parametersfixed<-c(parametersfixed, Tau=1)}

# Si Phi ou Delta sont indiqués en paramètres à ajuster, je les retire
parametersfit<-parametersfit[!names(parametersfit)=="Phi"]
parametersfit<-parametersfit[!names(parametersfit)=="Delta"]

# mais je les mets en paramètres fixes
if (is.na(parametersfixed["Phi"])) {parametersfixed<-c(parametersfixed, Phi=0)}
if (is.na(parametersfixed["Delta"])) {parametersfixed<-c(parametersfixed, Delta=0)}

	
if (progressbar) pb<-txtProgressBar(min=0, max=LDelta, style=3)

#parpre1<-parametersfit
parpre<-parametersfit

#FILLING MATRIX
for(j in 1:LDelta) {

XDelta<-Deltavalue[j]
for(i in 1:LPhi) {
  XPhi<-Phivalue[i]
  if (XDelta>=XPhi/2) {
    input[i,j]=NA
  } else {
  	parametersfixed["Delta"]<-XDelta
  	parametersfixed["Phi"]<-XPhi

#	assign("fixed", parametersfixed, envir=as.environment(.phenology.env))
    
    par<-parpre

    repeat {
    	    	
		resul<-optim(par, .Lnegbin, pt=list(data=data, fixed=parametersfixed, incertitude=method_incertitude, zerocounts=zero_counts), method="BFGS",control=list(trace=0, REPORT=1, maxit=500),hessian=FALSE)
		if (resul$convergence==0) break
		par<-resul$par
		# print("Convergence is not achieved. Optimization continues !")
	}
	
#	if (resul$value>378.717) {
#		print(resul$value)
#		print(parpre)
#		print(parametersfixed)
#	}
	
#  parpre<-resul$par
#  if (i==1) {parpre1<-parpre}
  input[i,j]<-resul$value
  }
#  Sys.sleep(0)
if (progressbar) setTxtProgressBar(pb, j)
}
}


Dv=as.vector(input)
pos=which.min(Dv)-1
j0=floor(pos/nrow(input))+1
i0=pos%%nrow(input)+1
print(paste("The minimum -Ln likelihood is ", input[i0, j0], sep=""))
print(paste("For Phi=",Phivalue[i0],sep=""))
print(paste("And Delta=",Deltavalue[j0],sep=""))

outputmap <- list(input=input, Phi=Phivalue, Delta=Deltavalue, Parametersfitted=names(parametersfit), 
Parametersfixed=parametersfixed, Data=names(data))

class(outputmap) <- "phenologymap"

return(outputmap)

}
}
