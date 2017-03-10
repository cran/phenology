#' par_init calculates initial set of parameters.
#' @title Calculate initial set of parameters.
#' @author Marc Girondot
#' @return The initial set of parameters
#' @param data Dataset generated with add_phenology()
#' @param fixed.parameters Set of fixed parameters
#' @param add.cofactors Names of cofactors that will be used (see fit_phenology)
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
#' Gratiot <- read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot <- add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg <- par_init(data_Gratiot, fixed.parameters=NULL)
#' # Run the optimisation
#' result_Gratiot <- fit_phenology(data=data_Gratiot, 
#' 		fitted.parameters=parg, fixed.parameters=NULL, trace=1)
#' data(result_Gratiot)
#' # Plot the phenology and get some stats
#' output<-plot(result_Gratiot)
#' }
#' @export


par_init <-
function(data=stop("A dataset must be provided"), 
         fixed.parameters=NULL, add.cofactors=NULL) {


if (is.null(fixed.parameters)) {fixed.parameters <- NA}

if (class(data)!="phenologydata") {
  stop("Data must be formated first using the function add_phenology().")
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

if (is.na(fixed.parameters["Min"]) & is.na(fixed.parameters["PMin"])) {

mB1<-0.1
names(mB1)<-paste("MinB_", names(data[serie]), sep="")

mE1<-0.1
names(mE1)<-paste("MinE_", names(data[serie]), sep="")

}

if (!is.na(fixed.parameters["Peak"])) {
	pkp<-fixed.parameters["Peak"]
} else {
	pkp<-od[which.max(nb)]
	pkp<-rnorm(1, pkp, 5)
	names(pkp)<-"Peak"
}

if (!is.na(fixed.parameters["Begin"])) {
	bgp<-fixed.parameters["Begin"]
} else {
	bgp<-od[1]+(pkp-od[1])/2
	names(bgp)<-"Begin"
}

if (!is.na(fixed.parameters["End"])) {
	edp<-fixed.parameters["End"]
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
if (is.na(fixed.parameters[c])) {par<-c(par, mx1)}
if (is.na(fixed.parameters["Min"]) & is.na(fixed.parameters["PMin"])) {
	c<-paste("MinB_", names(data[serie]), sep="")
	if ((is.na(fixed.parameters["MinB"]))*(is.na(fixed.parameters[c]))) {par<-c(par, mB1)}
	c<-paste("MinE_", names(data[serie]), sep="")
	if ((is.na(fixed.parameters["MinE"]))*(is.na(fixed.parameters[c]))) {par<-c(par, mE1)}
}

# fin de la boucle des séries

}

if ((is.na(fixed.parameters["Begin"])) && (is.na(fixed.parameters["Length"])) && (is.na(fixed.parameters["LengthB"]))) {par<-c(par, Begin=(min(bg)+mean(bg))/2)}
if ((is.na(fixed.parameters["Peak"]))) {par<-c(par, Peak=mean(pk))}
if ((is.na(fixed.parameters["End"])) && (is.na(fixed.parameters["Length"])) && (is.na(fixed.parameters["LengthE"]))) {par<-c(par, End=(max(ed)+mean(ed))/2)}
if (is.na(fixed.parameters["Flat"])) {par<-c(par, Flat=2)}

if (is.na(fixed.parameters["Theta"])) {
  par <- c(par, Theta=5)
}


par <- BE_to_LBLE(par)

# 19/3/2016
cof <- rep(0, length(add.cofactors))
names(cof) <- add.cofactors

cof <- cof[!(names(cof) %in% names(fixed.parameters))]

return(c(par, cof))
}
