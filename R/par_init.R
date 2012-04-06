#' par_init calculates initial set of parameters.
#' @title Calculate initial set of parameters.
#' @author Marc Girondot
#' @return The initial set of parameters
#' @param data dataset generated with add_format
#' @param parametersfixed Set of fixed parameters
#' @param help If TRUE, an help is displayed
#' @description This function is used to generate a first set of parameters
#' that is expected to be not to far from the final.
#' @export


par_init <-
function(data=dta, parametersfixed=NULL, help=FALSE) {

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

par<-NULL

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
	pk<-parametersfixed["Peak"]
} else {
	pk<-od[which.max(nb)]
	names(pk)<-"Peak"
}

if (!is.na(parametersfixed["Begin"])) {
	bg<-parametersfixed["Begin"]
} else {
	bg<-od[1]+(pk-od[1])/2
	names(bg)<-"Begin"
}

if (!is.na(parametersfixed["End"])) {
	bg<-parametersfixed["End"]
} else {
	ed<-od[length(od)]-(od[length(od)]-pk)/2
	names(ed)<-"End"
}

if ((pk<=bg)||(pk>=ed)) {
	pk<-bg+(ed-bg)/2
	names(pk)<-"Peak"
	}

#c(bg, pk, ed, Flat=2, mx1, mB1, mE1)

if ((serie==1) && (is.na(parametersfixed["Begin"])) && (is.na(parametersfixed["LengthB"]))) {par<-c(par, bg)}
if ((serie==1) && (is.na(parametersfixed["Peak"]))) {par<-c(par, pk)}
if ((serie==1) && (is.na(parametersfixed["End"])) && (is.na(parametersfixed["LengthE"]))) {par<-c(par, ed)}
if ((serie==1) && (is.na(parametersfixed["Flat"]))) {par<-c(par, Flat=2)}
c<-paste("Max_", names(data[serie]), sep="")
if (is.na(parametersfixed[c])) {par<-c(par, mx1)}
if (is.na(parametersfixed["Min"])) {
	c<-paste("MinB_", names(data[serie]), sep="")
	if ((is.na(parametersfixed["MinB"]))*(is.na(parametersfixed[c]))) {par<-c(par, mB1)}
	c<-paste("MinE_", names(data[serie]), sep="")
	if ((is.na(parametersfixed["MinE"]))*(is.na(parametersfixed[c]))) {par<-c(par, mE1)}
}

}

par<-c(par, Theta=5)

par<-BE_to_LBLE(par)

return(par)
}

}
