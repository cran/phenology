.daily_count <- function(d                  , 
                         xpar               , 
                         cofactors=NULL     , 
                         add.cofactors=NULL , 
                         print=FALSE        , 
                         zero=1E-9          ) {
  
  
  
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
  if (length(unique(index)) != 0) index <- as.character(unique(index)) else stop("Error in parameters\n", as.character(dput(xpar)))
  
  maxd <- 1+(1+max(d)%/%365)*365
  nn_m <- matrix(NA, ncol = maxd, nrow=length(index), dimnames = list(index, as.character((1:maxd)-1)))
  posd <- match(as.character(d), colnames(nn_m))
  dd <- (1:maxd)-1
  for (i in index) {
    
    idMinB <- (dd < xpar[paste0("Begin.", i)])
    idIncB <- (dd < xpar[paste0("PmoinsF.", i)]) & (!idMinB)
    idFlat <- (dd >= xpar[paste0("PmoinsF.", i)]) & (dd < xpar[paste0("PplusF.", i)])
    idIncE <- (dd >= xpar[paste0("PplusF.", i)]) & (dd < xpar[paste0("End.", i)])
    idMinE <- (dd >= xpar[paste0("End.", i)])
    
    nn_m[i, idMinB]  <- xpar[paste0("MinB.", i)]
    nn_m[i, idIncB]  <- ((1+cos(pi*(xpar[paste0("PmoinsF.", i)]-dd[idIncB])/xpar[paste0("PmoinsFB.", i)]))/2)*xpar[paste0("MaxMinB.", i)]+xpar[paste0("MinB.", i)]
    nn_m[i, idFlat]  <- xpar[paste0("Max.", i)]
    nn_m[i, idIncE]  <- ((1+cos(pi*(dd[idIncE]-(xpar[paste0("PplusF.", i)]))/xpar[paste0("EPplusF.", i)]))/2)*xpar[paste0("MaxMinE.", i)]+xpar[paste0("MinE.", i)]
    nn_m[i, idMinE]  <- xpar[paste0("MinE.", i)]
    ns <- ns1 <- ns2 <- ns3 <- rep(0, ncol(nn_m))
    if (xpar[paste0("sin.", i)]) {
      ns <- sin(2*pi*((dd+xpar[paste0("Delta.", i)])/xpar[paste0("Phi.", i)]))*(xpar[paste0("Alpha.", i)]+(xpar[paste0("Beta.", i)]*nn_m[i, ]^xpar[paste0("Tau.", i)]))
    }
    if (xpar[paste0("sin1.", i)]) {
      ns1 <- sin(2*pi*((dd+xpar[paste0("Delta1.", i)])/xpar[paste0("Phi1.", i)]))*(xpar[paste0("Alpha1.", i)]+(xpar[paste0("Beta1.", i)]*nn_m[i, ]^xpar[paste0("Tau1.", i)]))
    }
    if (xpar[paste0("sin2.", i)]) {
      ns2 <- sin(2*pi*((dd+xpar[paste0("Delta2.", i)])/xpar[paste0("Phi2.", i)]))*(xpar[paste0("Alpha2.", i)]+(xpar[paste0("Beta2.", i)]*nn_m[i, ]^xpar[paste0("Tau2.", i)]))
    }
    if (xpar[paste0("sin3.", i)]) {
      ns3 <- sin(2*pi*((dd+xpar[paste0("Delta3.", i)])/xpar[paste0("Phi3.", i)]))*(xpar[paste0("Alpha3.", i)]+(xpar[paste0("Beta3.", i)]*nn_m[i, ]^xpar[paste0("Tau3.", i)]))
    }
    nn_m[i, ] <- nn_m[i, ] + ns + ns1 + ns2 + ns3
  }
  
  # C'est quoi ça ??? 
  # C'est pour ne pas avoir de 0 mais pas sur xpar. Bizarre
  # xpar[paste0("PmoinsFB.", as.character(1:as.numeric(index)))] <- ifelse(xpar[paste0("PmoinsFB.", as.character(1:as.numeric(index)))] == 0, 
  #                                                            zero, 
  #                                                            xpar[paste0("PmoinsFB.", as.character(1:as.numeric(index)))])
  # xpar[paste0("EPplusF.", as.character(1:as.numeric(index)))] <- ifelse(xpar[paste0("EPplusF.", as.character(1:as.numeric(index)))] == 0, 
  #                                                           zero, 
  #                                                           xpar[paste0("EPplusF.", as.character(1:as.numeric(index)))])

  nn <- colSums(nn_m, dims=1)
  nn[is.na(nn)] <- zero
  nn[nn < zero] <- zero
  
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
  
  
  # je suis en en mode interactif, j'affiche le résultat
  if (print) {
    print(paste("Day ", d, "Number ", nn[posd]))
  }
  
  return(nn[posd])
}

