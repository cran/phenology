.daily_count <- function(d, xpar, 
                         cofactors=NULL, add.cofactors=NULL, 
                         print=FALSE, zero=1E-9) {
  

  
  # daily_count estimates nest number based on set of parameters.
  # @title Estimate expected counts based on set of parameters.
  # @author Marc Girondot
  # @param d Ordinal date (origin = 0)
  # @param xpar Set of fixed+fitted parameters
  # @param print If TRUE, the result is printed
  # @return The number of each day in d
  # @description Function estimates counts based on set of parameters.
  
  # Si c'est le modèle de Godley
  
  if (any(names(xpar) == "alpha", na.rm=TRUE)) {
    alpha <- xpar["alpha"]
    tp <- xpar["tp"]
    tf <- xpar["tf"]
    if (is.na(tf)) tf <- 0
    s1 <- xpar["s1"]
    s2 <- xpar["s2"]
    sr <- xpar["sr"]
    if (!is.na(sr)) {
      if (is.na(s1)) s1 <- sr
      if (is.na(s2)) s2 <- sr
    }
    
    return(ifelse(d<(tp-tf), alpha*exp(-((d-tp+tf)/s1)^2), 
                  ifelse(d>(tp+tf), alpha*exp(-((d-tp-tf)/s2)^2), 
                         alpha)))
  }
  
  # 10/5/2023
  index <- na.omit(suppressWarnings(as.numeric(gsub(".+\\.(\\d+).*", "\\1", names(xpar)))))
  if (length(index) != 0) index <- max(index) else index <- 0
  
  nn <- NULL
  
  if (!is.na(xpar["Begin"]))
  nn <- ifelse(d<xpar["Begin"], xpar["MinB"],
               ifelse(d<xpar["PmoinsF"], ((1+cos(pi*(xpar["PmoinsF"]-d)/xpar["PmoinsFB"]))/2)*xpar["MaxMinB"]+xpar["MinB"],
                      ifelse(d<xpar["PplusF"], xpar["Max"],
                             ifelse(d<xpar["End"], ((1+cos(pi*(d-(xpar["PplusF"]))/xpar["EPplusF"]))/2)*xpar["MaxMinE"]+xpar["MinE"],
                                    xpar["MinE"]
                             )
                      )
               )
  )
  
  if (index != 0)
    for (i in as.character(1:index)) {
      nn <- c(nn,  ifelse(d<xpar[paste0("Begin.", i)], xpar[paste0("MinB.", i)],
                     ifelse(d<xpar[paste0("PmoinsF.", i)], ((1+cos(pi*(xpar[paste0("PmoinsF.", i)]-d)/xpar[paste0("PmoinsFB.", i)]))/2)*xpar[paste0("MaxMinB.", i)]+xpar[paste0("MinB.", i)],
                            ifelse(d<xpar[paste0("PplusF.", i)], xpar[paste0("Max.", i)],
                                   ifelse(d<xpar[paste0("End.", i)], ((1+cos(pi*(d-(xpar[paste0("PplusF.", i)]))/xpar[paste0("EPplusF.", i)]))/2)*xpar[paste0("MaxMinE.", i)]+xpar[paste0("MinE.", i)],
                                          xpar[paste0("MinE.", i)]
                                   )
                            )
                     )
      )
      )
    }
  
  
  if (is.null(nn)) {
    print("No global model: Error, the parameters at the time of error are:")
    stop(dput(xpar))
  } else {
  if (any(is.na(nn))) {
    print("Global is NA: Error, the parameters at the time of error are:")
    stop(dput(xpar))
    #	assign("par_error", xpar, envir=as.environment(.phenology.env))
  }
  }
  
  nn <- matrix(nn, ncol = length(d), byrow = TRUE)
  nn <- colSums(nn, dims=1)
  
  
  if (xpar["sin"]) {
    if (xpar["Phi"] == 0) xpar["Phi"] <- 1E-9
    ns <- sin(2*pi*((d+xpar["Delta"])/xpar["Phi"]))*(xpar["Alpha"]+(xpar["Beta"]*nn^xpar["Tau"]))
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
    if (xpar["Phi1"] == 0) xpar["Phi1"] <- 1E-9
    ns1 <- sin(2*pi*((d+xpar["Delta1"])/xpar["Phi1"]))*(xpar["Alpha1"]+(xpar["Beta1"]*nn^xpar["Tau1"]))
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
    if (xpar["Phi2"] == 0) xpar["Phi2"] <- 1E-9
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
  
  nn[is.na(nn)] <- zero
  
  # Cofacteurs
  if ((!is.null(cofactors)) & (!is.null(add.cofactors))) {
    # Donne les paramètres cofacteurs
    xparcf <- xpar[(names(xpar) %in% add.cofactors) | (names(xpar) %in% paste0(add.cofactors, "multi"))]
    allxparcf <- rep(0, 2*length(add.cofactors))
    names(allxparcf) <- c(add.cofactors, paste0(add.cofactors, "multi"))
    xparcf <- modifyVector(val=xparcf, x=allxparcf)
    # cofactors$Date est une date
    # d est un nombre qui commence à 0
    # J'avais data$Date[i]
    effet1 <- rowSums(cofactors[cofactors$Date == d, add.cofactors, drop=FALSE] * xparcf[add.cofactors])
    effet2 <- rowSums(nn * cofactors[cofactors$Date == d, add.cofactors, drop=FALSE] * xparcf[paste0(add.cofactors, "multi")])
    nn <- nn + effet1 + effet2
  }

  nn[nn <= zero] <- zero
  nn[is.na(nn)] <- zero
  nn <- unname(nn)
  
  # je suis en en mode interactif, j'affiche le résultat
  if (print) {
    print(paste("Day ", d, "Number ", nn))
  }
  
  return(nn)
}

