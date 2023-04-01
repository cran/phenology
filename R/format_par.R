# .format_par calculates simplified set of parameters.
# @title The function ".format_par" is for internal use use only.
# @author Marc Girondot
# @return Return a modified set of parameters
# @param xpar Set of parameters
# @param serie Name of the series to be analyzed
# @description Calculate a simplified set of parameters.



.format_par <- function(xpar, serie) {
  
  #  xpar <- c(Min=12, Peak_Alpha=15, Peak_Beta=-16, Theta=16, Begin=15);serie="Alpha"
  #  getFromNamespace(".format_par", ns="phenology")(xpar, serie)
  
  # xpar <- c(Min=12, Peak_Alpha=15, Peak_Beta=-16, Theta_Alpha=16, Begin=15);serie="Alpha"
  # getFromNamespace(".format_par", ns="phenology")(xpar, serie)
  
  # xpar <- na.omit(xpar)
  
  
  #  xpar_courant <<- xpar
  #  save.image("courant.RData")
  
  nxparec <- strsplit(names(xpar), "_")
  # 1/11/2021 Le grepl me donne des résultats bizarres
  # gsub("([().])", "\\\\\\1", x[[2]]): je retire les (). et remplace par \\( et \\) et \\.
  ec <- sapply(nxparec, function(x) ifelse(length(x)>1, grepl(x[[2]], serie, fixed = TRUE), TRUE))
  # ec2 <- sapply(nxparec, function(x) ifelse(length(x)>1, substr(x[[2]], 1, nchar(serie))==serie, TRUE))
  # ec <- ec | ec2
  
  xparec <- xpar[ec]
  names(xparec) <- sapply(nxparec[ec], function(x) x[[1]])
  
  if (is.na(xparec["MinB"]) && is.na(xparec["PMinB"]) && is.na(xparec["Min"]) && is.na(xparec["PMin"])) {xparec["MinB"] <- 0}
  if (is.na(xparec["MinE"]) && is.na(xparec["PMinE"]) && is.na(xparec["Min"]) && is.na(xparec["PMin"])) {xparec["MinE"] <- 0}
  # 14/11/2018
  if (is.na(xparec["MinB.1"]) && is.na(xparec["PMinB.1"]) && is.na(xparec["Min.1"]) && is.na(xparec["PMin.1"])) {xparec["MinB.1"] <- 0}
  if (is.na(xparec["MinE.1"]) && is.na(xparec["PMinE.1"]) && is.na(xparec["Min.1"]) && is.na(xparec["PMin.1"])) {xparec["MinE.1"] <- 0}
  # 14/11/2018
  if (is.na(xparec["MinB.2"]) && is.na(xparec["PMinB.2"]) && is.na(xparec["Min.2"]) && is.na(xparec["PMin.2"])) {xparec["MinB.2"] <- 0}
  if (is.na(xparec["MinE.2"]) && is.na(xparec["PMinE.2"]) && is.na(xparec["Min.2"]) && is.na(xparec["PMin.2"])) {xparec["MinE.2"] <- 0}
  
  ec <- !is.na(match(names(xparec), 
                     c("MinB", "MinE", "Min", "PMin", "PMinE", "PMinB", "Peak", 
                       "Flat", "Flat.1", "Flat.2", 
                       "Begin", "End", "Max", "Theta", "Length", "LengthB", "LengthE", 
                       "Tau", "Tau1", "Tau2", 
                       "MinB.1", "MinE.1", "Min.1", "PMin.1", "PMinE.1", "PMinB.1", "Peak.1", 
                       "Begin.1", "End.1", "Max.1", "Length.1", "LengthB.1", "LengthE.1", 
                       "MinB.2", "MinE.2", "Min.2", "PMin.2", "PMinE.2", "PMinB.2", "Peak.2", 
                       "Begin.2", "End.2", "Max.2", "Length.2", "LengthB.2", "LengthE.2", 
                       "alpha", "tp", "s", "tf", "s1", "s2", "sr")))
  xparec[ec] <- abs(xparec[ec])
  
  xparec["Flat"] <- ifelse(is.na(xparec["Flat"]), 0, abs(xparec["Flat"]))
  xparec["Flat.1"] <- ifelse(is.na(xparec["Flat.1"]), 0, abs(xparec["Flat.1"]))
  xparec["Flat.2"] <- ifelse(is.na(xparec["Flat.2"]), 0, abs(xparec["Flat.2"]))
  
  if (!is.na(xparec["Length"])) {
    xparec["Begin"]=xparec["Peak"]-xparec["Length"]
    xparec["End"]=xparec["Peak"]+xparec["Length"]	
  }
  if (!is.na(xparec["LengthB"])) {
    xparec["Begin"]=xparec["Peak"]-xparec["LengthB"]
  }
  if (!is.na(xparec["LengthE"])) {
    xparec["End"]=xparec["Peak"]+xparec["LengthE"]	
  }
  
  if (!is.na(xparec["Length.1"])) {
    xparec["Begin.1"]=xparec["Peak.1"]-xparec["Length.1"]
    xparec["End.1"]=xparec["Peak.1"]+xparec["Length.1"]	
  }
  if (!is.na(xparec["LengthB.1"])) {
    xparec["Begin.1"]=xparec["Peak.1"]-xparec["LengthB.1"]
  }
  if (!is.na(xparec["LengthE.1"])) {
    xparec["End.1"]=xparec["Peak.1"]+xparec["LengthE.1"]	
  }
  
  if (!is.na(xparec["Length.2"])) {
    xparec["Begin.2"]=xparec["Peak.2"]-xparec["Length.2"]
    xparec["End.2"]=xparec["Peak.2"]+xparec["Length.2"]	
  }
  if (!is.na(xparec["LengthB.2"])) {
    xparec["Begin.2"]=xparec["Peak.2"]-xparec["LengthB.2"]
  }
  if (!is.na(xparec["LengthE.2"])) {
    xparec["End.2"]=xparec["Peak.2"]+xparec["LengthE.2"]	
  }
  
  
  
  if (!is.na(xparec["PMinE"])) {xparec["MinE"]<-xparec["Max"]*xparec["PMinE"]/100}
  if (!is.na(xparec["PMinB"])) {xparec["MinB"]<-xparec["Max"]*xparec["PMinB"]/100}

  if (!is.na(xparec["PMinE.1"])) {xparec["MinE.1"]<-xparec["Max.1"]*xparec["PMinE.1"]/100}
  if (!is.na(xparec["PMinB.1"])) {xparec["MinB.1"]<-xparec["Max.1"]*xparec["PMinB.1"]/100}

  if (!is.na(xparec["PMinE.2"])) {xparec["MinE.2"]<-xparec["Max.2"]*xparec["PMinE.2"]/100}
  if (!is.na(xparec["PMinB.2"])) {xparec["MinB.2"]<-xparec["Max.2"]*xparec["PMinB.2"]/100}
  
  if (!is.na(xparec["Min"])) {
    xparec["MinB"]<-xparec["Min"]
    xparec["MinE"]<-xparec["Min"]
  }
  
  if (!is.na(xparec["Min.1"])) {
    xparec["MinB.1"]<-xparec["Min.1"]
    xparec["MinE.1"]<-xparec["Min.1"]
  }
  
  if (!is.na(xparec["Min.2"])) {
    xparec["MinB.2"]<-xparec["Min.2"]
    xparec["MinE.2"]<-xparec["Min.2"]
  }
  
  if (!is.na(xparec["PMin"])) {
    xparec["MinB"]<-xparec["Max"]*xparec["PMin"]/100
    xparec["MinE"]<-xparec["Max"]*xparec["PMin"]/100	
  }
  
  if (!is.na(xparec["PMin.1"])) {
    xparec["MinB.1"]<-xparec["Max.1"]*xparec["PMin.1"]/100
    xparec["MinE.1"]<-xparec["Max.1"]*xparec["PMin.1"]/100	
  }
  
  if (!is.na(xparec["PMin.2"])) {
    xparec["MinB.2"]<-xparec["Max.2"]*xparec["PMin.2"]/100
    xparec["MinE.2"]<-xparec["Max.2"]*xparec["PMin.2"]/100	
  }
  
  xparec["sin"]<-(!is.na(xpar["Phi"]) && !is.na(xpar["Delta"]))
  if (xparec["sin"]) {
    if (is.na(xparec["Alpha"])) {xparec["Alpha"]=0}
    if (is.na(xparec["Beta"])) {xparec["Beta"]=0}
    if (is.na(xparec["Tau"])) {xparec["Tau"]=1}
  }
  
  xparec["sin1"]<-(!is.na(xpar["Phi1"]) && !is.na(xpar["Delta1"]))
  if (xparec["sin1"]) {
    if (is.na(xparec["Alpha1"])) {xparec["Alpha1"]=0}
    if (is.na(xparec["Beta1"])) {xparec["Beta1"]=0}
    if (is.na(xparec["Tau1"])) {xparec["Tau1"]=1}
  }
  
  xparec["sin2"]<-(!is.na(xpar["Phi2"]) && !is.na(xpar["Delta2"]))
  if (xparec["sin2"]) {
    if (is.na(xparec["Alpha2"])) {xparec["Alpha2"]=0}
    if (is.na(xparec["Beta2"])) {xparec["Beta2"]=0}
    if (is.na(xparec["Tau2"])) {xparec["Tau2"]=1}
  }
  
  
  xparec["PmoinsF"]<-xparec["Peak"]-(xparec["Flat"]/2)
  xparec["PplusF"]<-xparec["Peak"]+(xparec["Flat"]/2)
  
  xparec["PmoinsFB"]<-xparec["PmoinsF"]-xparec["Begin"]
  xparec["EPplusF"]<-xparec["End"]-xparec["PplusF"]
  
  xparec["MaxMinB"]<-xparec["Max"]-xparec["MinB"]
  xparec["MaxMinE"]<-xparec["Max"]-xparec["MinE"]
  
  # 1
  
  xparec["PmoinsF.1"]<-xparec["Peak.1"]-(xparec["Flat.1"]/2)
  xparec["PplusF.1"]<-xparec["Peak.1"]+(xparec["Flat.1"]/2)
  
  xparec["PmoinsFB.1"]<-xparec["PmoinsF.1"]-xparec["Begin.1"]
  xparec["EPplusF.1"]<-xparec["End.1"]-xparec["PplusF.1"]
  
  xparec["MaxMinB.1"]<-xparec["Max.1"]-xparec["MinB.1"]
  xparec["MaxMinE.1"]<-xparec["Max.1"]-xparec["MinE.1"]
  
  # 2
  
  xparec["PmoinsF.2"]<-xparec["Peak.2"]-(xparec["Flat.2"]/2)
  xparec["PplusF.2"]<-xparec["Peak.2"]+(xparec["Flat.2"]/2)
  
  xparec["PmoinsFB.2"]<-xparec["PmoinsF.2"]-xparec["Begin.2"]
  xparec["EPplusF.2"]<-xparec["End.2"]-xparec["PplusF.2"]
  
  xparec["MaxMinB.2"]<-xparec["Max.2"]-xparec["MinB.2"]
  xparec["MaxMinE.2"]<-xparec["Max.2"]-xparec["MinE.2"]
  

  
  return(xparec)
  
}
