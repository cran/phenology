
.daily_count <- function(d, xpar, print=FALSE, zero=1E-9) {
  
  # daily_count estimates nest number based on set of parameters.
  # @title Estimate expected counts based on set of parameters.
  # @author Marc Girondot
  # @param d Ordinal date
  # @param xpar Set of fixed+fitted parameters
  # @param print If TRUE, the result is displayed
  # @param help If TRUE, an help is displayed
  # @return The number of each day in d
  # @description Function estimates counts based on set of parameters.
  
  
  nn <- ifelse(d<xpar["Begin"], xpar["MinB"],
               ifelse(d<xpar["PmoinsF"], ((1+cos(pi*(xpar["PmoinsF"]-d)/xpar["PmoinsFB"]))/2)*xpar["MaxMinB"]+xpar["MinB"],
                      ifelse(d<xpar["PplusF"], xpar["Max"],
                             ifelse(d<xpar["End"], ((1+cos(pi*(d-(xpar["PplusF"]))/xpar["EPplusF"]))/2)*xpar["MaxMinE"]+xpar["MinE"],
                                    xpar["MinE"]
                             )
                      )
               )
  )
  
  nn.1 <- ifelse(d<xpar["Begin.1"], xpar["MinB.1"],
               ifelse(d<xpar["PmoinsF.1"], ((1+cos(pi*(xpar["PmoinsF.1"]-d)/xpar["PmoinsFB.1"]))/2)*xpar["MaxMinB.1"]+xpar["MinB.1"],
                      ifelse(d<xpar["PplusF.1"], xpar["Max.1"],
                             ifelse(d<xpar["End.1"], ((1+cos(pi*(d-(xpar["PplusF.1"]))/xpar["EPplusF.1"]))/2)*xpar["MaxMinE.1"]+xpar["MinE.1"],
                                    xpar["MinE.1"]
                             )
                      )
               )
  )
  
  nn.2 <- ifelse(d<xpar["Begin.2"], xpar["MinB.2"],
               ifelse(d<xpar["PmoinsF.2"], ((1+cos(pi*(xpar["PmoinsF.2"]-d)/xpar["PmoinsFB.2"]))/2)*xpar["MaxMinB.2"]+xpar["MinB.2"],
                      ifelse(d<xpar["PplusF.2"], xpar["Max.2"],
                             ifelse(d<xpar["End.2"], ((1+cos(pi*(d-(xpar["PplusF.2"]))/xpar["EPplusF.2"]))/2)*xpar["MaxMinE.2"]+xpar["MinE.2"],
                                    xpar["MinE.2"]
                             )
                      )
               )
  )
  
  nn <- ifelse(is.na(nn), 0, nn)
  nn.1 <- ifelse(is.na(nn.1), 0, nn.1)
  nn.2 <- ifelse(is.na(nn.2), 0, nn.2)
  nn <- nn + nn.1 + nn.2
  
  if (any(is.na(nn))) {
    print("Global: Error, the parameters at the time of error are:")
    stop(dput(xpar))
    #	assign("par_error", xpar, envir=as.environment(.phenology.env))
  }
  
  
  if (xpar["sin"]) {
    ns<-sin(2*pi*((d+xpar["Delta"])/xpar["Phi"]))*(xpar["Alpha"]+(xpar["Beta"]*nn^xpar["Tau"]))
    if (any(is.na(ns))) {
      print(d)
      print("Sin: Error, the parameters at the time of error are:")
      stop(dput(xpar))
      #	assign("par_error", xpar, envir=as.environment(.phenology.env))
    }
  } else {
    ns <- 0
  }
  
  if (xpar["sin1"]) {
    ns1<-sin(2*pi*((d+xpar["Delta1"])/xpar["Phi1"]))*(xpar["Alpha1"]+(xpar["Beta1"]*nn^xpar["Tau1"]))
    if (any(is.na(ns1))) {
      print(d)
      print("Sin 1: Error, the parameters at the time of error are:")
      stop(dput(xpar))
      #	assign("par_error", xpar, envir=as.environment(.phenology.env))
    }
  } else {
    ns1 <- 0
  }
  
  
  if (xpar["sin2"]) {
    ns2<-sin(2*pi*((d+xpar["Delta2"])/xpar["Phi2"]))*(xpar["Alpha2"]+(xpar["Beta2"]*nn^xpar["Tau2"]))
    if (any(is.na(ns2))) {
      print(d)
      print("Sin 2: Error, the parameters at the time of error are:")
      stop(dput(xpar))
      #	assign("par_error", xpar, envir=as.environment(.phenology.env))
    }
  } else {
    ns2 <- 0
  }
  
  nn <- nn+ns+ns1+ns2
  nn[is.na(nn)] <- zero
  
  nn <- ifelse((nn <= zero) & (d<xpar["Begin"]), xpar["MinB"],
               ifelse((nn <= zero) & (d>xpar["Begin"]), xpar["MinE"],
                      ifelse(nn <= zero,(xpar["MinB"]+xpar["MinE"])/2,
                             nn
                      )
               )
  )
  
  nn[nn <= zero] <- zero
  nn[is.na(nn)] <- zero
  nn <- unname(nn)
  
  # je suis en en mode interactif, j'affiche le résultat
  if (print) {
    print(paste("Day ", d, "Number ", nn))
  }
  
  return(nn)
}

