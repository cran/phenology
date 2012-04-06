#' .format_par calculates simplified set of parameters.
#' @title The function ".format_par" is for internal use use only.
#' @author Marc Girondot
#' @return Return a modified set of parameters
#' @param xpar Set of parameters
#' @param serie Standard deviation value to be added
#' @param help If TRUE, an help is displayed
#' @description Calculate a simplified set of parameters.



.format_par <-
function(xpar, serie, help=FALSE) {
if(help) {
	cat("This function is for internal use use only.\n")
} else {


xparec<-NULL
# je ne stocke que les paramètres associés à la série en cours
	for(ii in 1:length(xpar)) {
		ec<-TRUE
		parec<-xpar[ii]
		if (substr(names(parec), 1, 4)=="Max_" && names(parec)!=paste("Max_", serie, sep="")) {ec<-FALSE}
		if (substr(names(parec), 1, 4)=="Min_" && names(parec)!=paste("Min_", serie, sep="")) {ec<-FALSE}
		if (substr(names(parec), 1, 5)=="MinE_" && names(parec)!=paste("MinE_", serie, sep="")) {ec<-FALSE}
		if (substr(names(parec), 1, 5)=="MinB_" && names(parec)!=paste("MinB_", serie, sep="")) {ec<-FALSE}
		if (substr(names(parec), 1, 5)=="Peak_" && names(parec)!=paste("Peak_", serie, sep="")) {ec<-FALSE}

		names(parec)<-strsplit(names(parec), "_")[[1]][1]

		if(ec) {xparec<-c(xparec, parec)}
		
	}

if (is.na(xparec["MinB"]) && is.na(xparec["PMinB"]) && is.na(xparec["Min"]) && is.na(xparec["PMin"])) {xparec["MinB"]<-0}
if (is.na(xparec["MinE"]) && is.na(xparec["PMinE"]) && is.na(xparec["Min"]) && is.na(xparec["PMin"])) {xparec["MinE"]<-0}

xparec["Flat"]<-ifelse(is.na(xparec["Flat"]), 0, abs(xparec["Flat"]))
xparec["Theta"]<-abs(xparec["Theta"])

if (!is.na(xparec["MinB"])) {xparec["MinB"]=abs(xparec["MinB"])}
if (!is.na(xparec["MinE"])) {xparec["MinE"]=abs(xparec["MinE"])}
if (!is.na(xparec["Min"])) {xparec["Min"]=abs(xparec["Min"])}
if (!is.na(xparec["Pmin"])) {xparec["PMin"]=abs(xparec["PMin"])}
if (!is.na(xparec["PminB"])) {xparec["PMinB"]=abs(xparec["PMinB"])}
if (!is.na(xparec["PminE"])) {xparec["PMinE"]=abs(xparec["PMinE"])}
if (!is.na(xparec["Max"])) {xparec["Max"]=abs(xparec["Max"])}
if (!is.na(xparec["Length"])) {
	xparec["Length"]=abs(xparec["Length"])
	xparec["Begin"]=xparec["Peak"]-xparec["Length"]
	xparec["End"]=xparec["Peak"]+xparec["Length"]	
	}
if (!is.na(xparec["LengthB"])) {
	xparec["LengthB"]=abs(xparec["LengthB"])
	xparec["Begin"]=xparec["Peak"]-xparec["LengthB"]
	}
if (!is.na(xparec["LengthE"])) {
	xparec["LengthE"]=abs(xparec["LengthE"])
	xparec["End"]=xparec["Peak"]+xparec["LengthE"]	
	}



if (!is.na(xparec["PMinE"])) {xparec["MinE"]<-xparec["Max"]*xparec["PMinE"]/100}
if (!is.na(xparec["PMinB"])) {xparec["MinB"]<-xparec["Max"]*xparec["PMinB"]/100}

if (!is.na(xparec["Min"])) {
	xparec["MinB"]<-xparec["Min"]
	xparec["MinE"]<-xparec["Min"]
	}
	
if (!is.na(xparec["PMin"])) {
	xparec["MinB"]<-xparec["Max"]*xparec["PMin"]/100
	xparec["MinE"]<-xparec["Max"]*xparec["PMin"]/100	
	}

xparec["sin"]<-(!is.na(xpar["Phi"]) && !is.na(xpar["Delta"]))
if (xparec["sin"]) {
	if (is.na(xparec["Alpha"])) {xparec["Alpha"]=0}
	if (is.na(xparec["Beta"])) {xparec["Beta"]=0}
	if (is.na(xparec["Tau"])) {xparec["Tau"]=1}
	xparec["Tau"]=abs(xparec["Tau"])
}

xparec["sin1"]<-(!is.na(xpar["Phi1"]) && !is.na(xpar["Delta1"]))
if (xparec["sin1"]) {
	if (is.na(xparec["Alpha1"])) {xparec["Alpha1"]=0}
	if (is.na(xparec["Beta1"])) {xparec["Beta1"]=0}
	if (is.na(xparec["Tau1"])) {xparec["Tau1"]=1}
	xparec["Tau1"]=abs(xparec["Tau1"])
}

xparec["sin2"]<-(!is.na(xpar["Phi2"]) && !is.na(xpar["Delta2"]))
if (xparec["sin2"]) {
	if (is.na(xparec["Alpha2"])) {xparec["Alpha2"]=0}
	if (is.na(xparec["Beta2"])) {xparec["Beta2"]=0}
	if (is.na(xparec["Tau2"])) {xparec["Tau2"]=1}
	xparec["Tau2"]=abs(xparec["Tau2"])
}


xparec["PmoinsF"]<-xparec["Peak"]-(xparec["Flat"]/2)
xparec["PplusF"]<-xparec["Peak"]+(xparec["Flat"]/2)

xparec["PmoinsFB"]<-xparec["PmoinsF"]-xparec["Begin"]
xparec["EPplusF"]<-xparec["End"]-xparec["PplusF"]

xparec["MaxMinB"]<-xparec["Max"]-xparec["MinB"]
xparec["MaxMinE"]<-xparec["Max"]-xparec["MinE"]

return(xparec)

}
}
