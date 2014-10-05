#' par_init calculates initial set of parameters.
#' @title Calculate initial set of parameters.
#' @author Marc Girondot
#' @return The initial set of parameters
#' @param data Dataset generated with add_phenology()
#' @param parametersfixed Set of fixed parameters
#' @param help If TRUE, an help is displayed
#' @description This function is used to generate a first set of parameters
#' that is expected to be not to far from the final.\cr
#' The parameters can be:\cr
#' Min, MinE, MinB, PMin, PMinB, PMinE\cr
#' Max\cr
#' Begin, Peak, Flat, End\cr
#' theta\cr
#' Alpha, Beta, tau, Phi, Delta\cr
#' Alpha1, Beta1, tau1, Phi1, Delta1\cr
#' Alpha2, Beta2, tau2, Phi2, Delta2\cr
#' Alpha3, Beta3, tau3, Phi3, Delta3\cr
#' The parameters Max, Min, MinE, MinB, and Peak can be followed with _ and the name of the rookery.
#' @examples
#' \dontrun{
#' library(phenology)
#' # Read a file with data
#' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' result_Gratiot<-fit_phenology(data=data_Gratiot, 
#' 		parametersfit=parg, parametersfixed=NULL, trace=1)
#' data(result_Gratiot)
#' # Plot the phenology and get some stats
#' output<-plot(result_Gratiot)
#' }
#' @export


par_init <-
function(data=stop("A dataset must be provided"), parametersfixed=NA, help=FALSE) {

if(help) {
	cat("This function is used to generate a first set of parameters\n")
	cat("that is expected to be not to far from the final\n")
	cat("parameters after fitting.\n")
	cat("The general syntax is:\n")
	cat("par<-par_init(data=dataset, parametersfixed=parfixed)\n")
	cat("with parfixed being the list of parameters that should\n")
	cat("not be fitted. The parameter names are:\n")
	cat("Begin Peak End Flat Max Min MinE MinB Theta\n")
	cat("For example: parfixed<-c(Min=0, Flat=0)\n")
	cat("means that both Min and Flat parameters must be\n")
	cat("forced to 0.\n")
	cat("To use the B and E formulation, you\n")
	cat("must convert this set of parameters using:\n")
	cat("par<-LBLE_to_BE(par)\n")
	cat("Use the function add_SE to add standard error to fixed parameters.\n")
} else {

if (is.null(parametersfixed)) {parametersfixed<-NA}

if (class(data)!="phenologydata") {
  cat("Data must be formated first using the function add_format().\n")
  return()
}


par<-NULL

bg <- NULL
pk <- NULL
ed <- NULL

for(serie in 1:length(data)) {

od<-data[[serie]]$ordinal
nb<-as.numeric(data[[serie]]$nombre)

mx1<-max(nb, na.rm=TRUE)/2
names(mx1)<-paste("Max_", names(data[serie]), sep="")

if (is.na(parametersfixed["Min"])) {

mB1<-0.5
names(mB1)<-paste("MinB_", names(data[serie]), sep="")

mE1<-0.5
names(mE1)<-paste("MinE_", names(data[serie]), sep="")

}

if (!is.na(parametersfixed["Peak"])) {
	pkp<-parametersfixed["Peak"]
} else {
	pkp<-od[which.max(nb)]
	pkp<-rnorm(1, pkp, 5)
	names(pkp)<-"Peak"
}

if (!is.na(parametersfixed["Begin"])) {
	bgp<-parametersfixed["Begin"]
} else {
	bgp<-od[1]+(pkp-od[1])/2
	names(bgp)<-"Begin"
}

if (!is.na(parametersfixed["End"])) {
	edp<-parametersfixed["End"]
} else {
	edp<-od[length(od)]-(od[length(od)]-pkp)/2
	names(edp)<-"End"
}

if ((pkp<=bgp)||(pkp>=edp)) {
	pkp<-bgp+(edp-bgp)/2	
	names(pkp)<-"Peak"
	}
	
bg <- c(bg, bgp)
ed <- c(ed, edp)
pk <- c(pk, pkp)

#c(bg, pk, ed, Flat=2, mx1, mB1, mE1)

c<-paste("Max_", names(data[serie]), sep="")
if (is.na(parametersfixed[c])) {par<-c(par, mx1)}
if (is.na(parametersfixed["Min"])) {
	c<-paste("MinB_", names(data[serie]), sep="")
	if ((is.na(parametersfixed["MinB"]))*(is.na(parametersfixed[c]))) {par<-c(par, mB1)}
	c<-paste("MinE_", names(data[serie]), sep="")
	if ((is.na(parametersfixed["MinE"]))*(is.na(parametersfixed[c]))) {par<-c(par, mE1)}
}

# fin de la boucle des séries

}

if ((is.na(parametersfixed["Begin"])) && (is.na(parametersfixed["Length"])) && (is.na(parametersfixed["LengthB"]))) {par<-c(par, Begin=(min(bg)+mean(bg))/2)}
if ((is.na(parametersfixed["Peak"]))) {par<-c(par, Peak=mean(pk))}
if ((is.na(parametersfixed["End"])) && (is.na(parametersfixed["Length"])) && (is.na(parametersfixed["LengthE"]))) {par<-c(par, End=(max(ed)+mean(ed))/2)}
if ((is.na(parametersfixed["Flat"]))) {par<-c(par, Flat=2)}


par<-c(par, Theta=5)

par<-BE_to_LBLE(par)

return(par)
}

}
