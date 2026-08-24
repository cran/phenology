# .format_par calculates simplified set of parameters.
# @title The function ".format_par" is for internal use use only.
# @author Marc Girondot
# @return Return a modified set of parameters
# @param xpar Set of parameters
# @param serie Name of the series to be analyzed
# @description Calculate a simplified set of parameters.

# xpar <- c('Max.1' = 29.80724298715009, 
#   'MinB.1' = 5,
#   'MinE.1' = 5,
#   'LengthB.1' = 48.874255673500919, 
#   'Peak.1' = 124.88392063302938, 
#   'LengthE.1' = 179.86385332359652, 
#   'Max.2' = 5.6693380329255865,
#   'MinB.2' = 5,
#   'MinE.2' = 5,
#   'LengthB.2' = 153.6648921722695, 
#   'Peak.2' = 345.96111524736443, 
#   'LengthE.2' = 0.048129714036260852, 
#   'Theta' = 10.584509187364276)

.format_par <- function(xpar, 
                        serie, 
                        model_before=NULL, 
                        season=NULL      ) {
  
  # if (!is.null(model_before)) eval(parse(text=model_before), envir= environment())
  # model_before <- "Peak.1=Peak.3; Max.1=Max.2"
  # model_before = "Length.3 = Length.1; Length.4= Length.2; Length.5 = Length.1"
  
  if (!is.null(model_before)) {
    model_before <- gsub(" ", "", model_before)
    model_before <- strsplit(model_before, ";")[[1]]
    for (i in 1:length(model_before)) {
      model_before[i] <- paste0("xpar['", gsub("=", "']=xpar['", model_before[i]), "']")
    }
    model_before <- paste0(model_before, collapse = ";")
    eval(parse(text=model_before), envir= environment())
  }
  
  #  xpar <- c(Min=12, Peak_Alpha=15, Peak_Beta=-16, Theta=16, Begin=15);serie="Alpha"
  #  getFromNamespace(".format_par", ns="phenology")(xpar, serie)
  
  # xpar <- c(Min=12, Peak_essai=15, Peak_Beta=-16, Theta_essai=16, Begin=15);serie="essai"
  # getFromNamespace(".format_par", ns="phenology")(xpar, serie)
  
  # xpar <- c(Min.2=12, Peak.1_essai=15, Peak.2_essai=-16, Theta_essai=16, Begin=15);serie="essai"
  
  
  # xpar <- na.omit(xpar)
  
  
  #  xpar_courant <<- xpar
  #  save.image("courant.RData")
  
  nxparec <- strsplit(names(xpar), "_")
  # Dans ec j'ai TRUE si je garde la série sur la base soit de son nom, soit pas de nom
  ec <- sapply(nxparec, function(x) ifelse(length(x)>1, grepl(x[[2]], serie, fixed = TRUE), TRUE))
  
  # Je prends ceux de la série en cours
  xparec <- xpar[ec]
  names(xparec) <- sapply(nxparec[ec], function(x) x[[1]])
  # J'ai deux fois le même paramètre après avoir retiré le nom
  if (length(unique(names(xparec))) != length(names(xparec))) stop("At least two series have similar names that can be confounded.")
  # 10/5/2023. Si j'ai un .nombre, j'ai un index
  index <- na.omit(suppressWarnings(as.numeric(gsub(".+\\.(\\d+).*", "\\1", names(xparec)))))
  if (length(index) != 0) index <- max(index) else index <- 0
  
  # A ce moment dans xparec je n'ai plus le nom de la série
  # je mets un index à tous
  if (index == 0) {
    names(xparec) <- paste0(names(xparec), ".1")
    index <- 1
  }
  
  # Si j'ai un paramètre sans .x, je le duplique index fois
  xparec_ec <- NULL
  for (nec in names(xparec)) {
    if (!grepl("\\.", nec)) {
      for (i in paste0(".", as.character(1:index))) {
        if (is.na(xparec[paste0(nec, i)])) {
          xparec_ec2 <- xparec[nec]
          names(xparec_ec2) <- paste0(names(xparec_ec2), i)
          xparec_ec <- c(xparec_ec, xparec_ec2)
        }
      }
    } else {
      # il faut aussi le dupliquer s'il manque des index
      xparec_ec <- c(xparec_ec, xparec[nec])
    }
  }
  xparec <- xparec_ec
  
  # les différents paramètres sont unique(gsub("\\.d+", "", names(xparec)))
  if (index >1) {
    xparec_ec <- NULL
    for (nec in unique(gsub("\\.[0-9]+", "", names(xparec)))) {
      if (length(sum(grepl(paste0("^", nec, "$"), gsub("\\.[0-9]+", "", names(xparec))))) < index) {
        if (length(sum(grepl(paste0("^", nec, "$"), gsub("\\.[0-9]+", "", names(xparec))))) > 1) {
          stop(paste0("I dn't know which ", nec, "to choose."))
        } else {
          for (i in paste0(".", as.character(1:index))) {
            if (is.na(xparec[paste0(nec, i)])) {
              xparec_ec2 <- xparec[names(xparec)[grepl(paste0("^", nec, "."), names(xparec))]]
              names(xparec_ec2) <- paste0(nec, i)
              xparec_ec <- c(xparec_ec, xparec_ec2)
            }
          }
        }
      }
    }
    xparec <- c(xparec, xparec_ec)
  }
  
  xparec_ec <- NULL
  # Maintenant je gère les périodes
  for (i in paste0(".", as.character(1:index))) {
    xparec_ec[paste0("sin", i)] <- (any(!is.na(xparec[grepl(paste0("^Phi", i, "$"), names(xparec))]))) & (any(!is.na(xparec[grepl(paste0("^Delta", i, "$"), names(xparec))])))
    if (xparec_ec[paste0("sin", i)]) {
      if (is.na(xparec[paste0("Alpha", i)])) {xparec_ec[paste0("Alpha", i)]=0}
      if (is.na(xparec[paste0("Beta", i)])) {xparec_ec[paste0("Beta", i)]=0}
      if (is.na(xparec[paste0("Tau", i)])) {xparec_ec[paste0("Tau", i)]=1}
    }
    xparec_ec[paste0("sin1", i)] <- (any(!is.na(xparec[grepl(paste0("^Phi1", i, "$"), names(xparec))]))) & (any(!is.na(xparec[grepl(paste0("^Delta1", i, "$"), names(xparec))])))
    if (xparec_ec[paste0("sin1", i)]) {
      if (is.na(xparec[paste0("Alpha1", i)])) {xparec_ec[paste0("Alpha1", i)]=0}
      if (is.na(xparec[paste0("Beta1", i)])) {xparec_ec[paste0("Beta1", i)]=0}
      if (is.na(xparec[paste0("Tau1", i)])) {xparec_ec[paste0("Tau1", i)]=1}
    }
    xparec_ec[paste0("sin2", i)] <- (any(!is.na(xparec[grepl(paste0("^Phi2", i, "$"), names(xparec))]))) & (any(!is.na(xparec[grepl(paste0("^Delta2", i, "$"), names(xparec))])))
    if (xparec_ec[paste0("sin2", i)]) {
      if (is.na(xparec[paste0("Alpha2", i)])) {xparec_ec[paste0("Alpha2", i)]=0}
      if (is.na(xparec[paste0("Beta2", i)])) {xparec_ec[paste0("Beta2", i)]=0}
      if (is.na(xparec[paste0("Tau2", i)])) {xparec_ec[paste0("Tau2", i)]=1}
    }
    xparec_ec[paste0("sin3", i)] <- (any(!is.na(xparec[grepl(paste0("^Phi3", i, "$"), names(xparec))]))) & (any(!is.na(xparec[grepl(paste0("^Delta3", i, "$"), names(xparec))])))
    if (xparec_ec[paste0("sin3", i)]) {
      if (is.na(xparec[paste0("Alpha3", i)])) {xparec_ec[paste0("Alpha3", i)]=0}
      if (is.na(xparec[paste0("Beta3", i)])) {xparec_ec[paste0("Beta3", i)]=0}
      if (is.na(xparec[paste0("Tau3", i)])) {xparec_ec[paste0("Tau3", i)]=1}
    }
    if (is.na(xparec[paste0("MinB", i)]) && is.na(xparec[paste0("PMinB", i)]) && is.na(xparec[paste0("Min", i)]) && is.na(xparec[paste0("PMin", i)])) {xparec_ec[paste0("MinB", i)] <- 0}
    if (is.na(xparec[paste0("MinE", i)]) && is.na(xparec[paste0("PMinE", i)]) && is.na(xparec[paste0("Min", i)]) && is.na(xparec[paste0("PMin", i)])) {xparec_ec[paste0("MinE", i)] <- 0}
    if (is.na(xparec[paste0("Flat", i)])) {xparec[paste0("Flat", i)] <- 0}
  }
  
  xparec <- c(xparec, xparec_ec)
  
  for (i in as.character(1:index)) {
    if (is.na(xparec[paste0("MinB.", i)]) && is.na(xparec[paste0("PMinB.", i)]) && is.na(xparec[paste0("Min.", i)]) && is.na(xparec[paste0("PMin.", i)])) {xparec[paste0("MinB.", i)] <- 0}
    if (is.na(xparec[paste0("MinE.", i)]) && is.na(xparec[paste0("PMinE.", i)]) && is.na(xparec[paste0("Min.", i)]) && is.na(xparec[paste0("PMin.", i)])) {xparec[paste0("MinE.", i)] <- 0}
    if (is.na(xparec[paste0("Flat.", i)])) {xparec[paste0("Flat.", i)] <- 0}
  }
  
  # Ca ne va pas
  
  ec <- grepl("^Min|^PMin|^Peak|^Flat|^Begin|^End|^Max|^Theta|^Length|^Tau|^alpha|^tp|^s|^tf|^s1|^s2|^sr", names(xparec))
  
  # ec <- !is.na(match(names(xparec), 
  #                    c("MinB", "MinE", "Min", "PMin", "PMinE", "PMinB", "Peak", 
  #                      "Flat", 
  #                      "Begin", "End", "Max", "Theta", "Length", "LengthB", "LengthE", 
  #                      "Tau", "Tau1", "Tau2", 
  #                      "alpha", "tp", "s", "tf", "s1", "s2", "sr")))
  xparec[ec] <- abs(xparec[ec])
  
  
  # if (index != 0)
  #   # for (i in as.character(1:index)) {
  #   #   if (is.na(xparec[paste0("Flat.", i)])) xparec[paste0("Flat.", i)] <- xparec["Flat"]
  #   #   # xparec[paste0("Flat.", i)] <- ifelse(is.na(xparec[paste0("Flat.", i)]), 0, abs(xparec[paste0("Flat.", i)]))
  #   # }
  
  
  # Je considère que Peak existe
  for (i in as.character(1:index)) {
    if (!is.na(xparec[paste0("Length.", i)]))  {
      xparec[paste0("LengthB.", i)] <- xparec[paste0("Length.", i)]
      xparec[paste0("LengthE.", i)] <- xparec[paste0("Length.", i)]
    }
    if (is.na(xparec[paste0("LengthB.", i)])) {
      xparec[paste0("LengthB.", i)] <- xparec[paste0("Peak.", i)] - xparec[paste0("Begin", i)]
    }
    if (is.na(xparec[paste0("LengthE.", i)])) {
      xparec[paste0("LengthE.", i)] <- xparec[paste0("End", i)] - xparec[paste0("Peak.", i)]
    }
    if (is.na(xparec[paste0("Begin.", i)])) {
      xparec[paste0("Begin.", i)] <- xparec[paste0("Peak.", i)] - xparec[paste0("LengthB.", i)]
    }
    if (is.na(xparec[paste0("End.", i)])) {
      xparec[paste0("End.", i)] <- xparec[paste0("Peak.", i)] + xparec[paste0("LengthE.", i)]
    }
    if (!is.na(xparec[paste0("PMinE.", i)])) {xparec[paste0("MinE.", i)]<-xparec[paste0("Max.", i)]*xparec[paste0("PMinE.", i)]/100}
    if (!is.na(xparec[paste0("PMinB.", i)])) {xparec[paste0("MinB.", i)]<-xparec[paste0("Max.", i)]*xparec[paste0("PMinB.", i)]/100}
    if (!is.na(xparec[paste0("PMin.", i)])) {
      xparec[paste0("MinB.", i)]<-xparec[paste0("Max.", i)]*xparec[paste0("PMin.", i)]/100
      xparec[paste0("MinE.", i)]<-xparec[paste0("Max.", i)]*xparec[paste0("PMin.", i)]/100
    }
    if (!is.na(xparec[paste0("Min.", i)])) {
      xparec[paste0("MinB.", i)] <- xparec[paste0("Min.", i)]
      xparec[paste0("MinE.", i)] <- xparec[paste0("Min.", i)]
    }
    
    xparec[paste0("PmoinsF.", i)] <- xparec[paste0("Peak.", i)]-(xparec[paste0("Flat.", i)]/2)
    xparec[paste0("PplusF.", i)] <- xparec[paste0("Peak.", i)]+(xparec[paste0("Flat.", i)]/2)
    
    xparec[paste0("PmoinsFB.", i)] <- xparec[paste0("PmoinsF.", i)]-xparec[paste0("Begin.", i)]
    xparec[paste0("EPplusF.", i)] <- xparec[paste0("End.", i)]-xparec[paste0("PplusF.", i)]
    
    xparec[paste0("MaxMinB.", i)] <- xparec[paste0("Max.", i)]-xparec[paste0("MinB.", i)]
    xparec[paste0("MaxMinE.", i)] <- xparec[paste0("Max.", i)]-xparec[paste0("MinE.", i)]
  }
  
  
  # Je ne garde qu'un seul Theta
  th <- xparec["Theta.1"]
  xparec <- xparec[!grepl("^Theta", names(xparec))]
  xparec["Theta"] <- unname(th)
  
  
  # JE ne suis pas sûr que c'est là qu'il faut le mettre
  if (!is.null(season)) {
    xparec_ec <- NULL
    for (i in as.character(season)) {
      xparec_i <- xparec[grepl(paste0("\\.", i,"$"), names(xparec))]
      xparec_ec <- c(xparec_ec, xparec_i)
    }
    xparec <- c(xparec_ec, xparec["Theta"])
  }
  
  return(xparec)
  
}
