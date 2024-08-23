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

.format_par <- function(xpar, serie, model_before=NULL, season=NULL) {
  
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
  
  # Je prends ceux de la série en cours
  xparec <- xpar[ec]
  names(xparec) <- sapply(nxparec[ec], function(x) x[[1]])
  
  if (length(unique(names(xparec))) != length(names(xparec))) stop("At least two series have similar names that can be confound")
  # 10/5/2023
  index <- na.omit(suppressWarnings(as.numeric(gsub(".+\\.(\\d+).*", "\\1", names(xparec)))))
  if (length(index) != 0) index <- max(index) else index <- 0
  
  # Je garde que les paramètres avec des . et je les crée si ils n'existent pas
  if (index != 0) {
    # Si j'ai plus d'un index
    # Les nouveaux paramètres que je crée
    xparec_ec <- NULL
    # Je les passe en revue un par un
    for (i in 1:length(xparec)) {
      # Si je n'ai pas de point, sauf Theta
      if ((!grepl("\\.", names(xparec[i]))) & (names(xparec[i]) != "Theta")) {
        # Je prends le paramètre sans point
        xx <- xparec[i]
        for (j in 1:index) {
          xxp <- xparec[i]
          names(xxp) <- paste0(names(xparec[i]), ".", as.character(j))
          xparec_ec <- c(xparec_ec, xxp)
        }
      } else {
        xparec_ec <- c(xparec_ec, xparec[i])
      }
    }
    xparec <- xparec_ec
    for (i in as.character(1:index)) {
      if (is.na(xparec[paste0("MinB.", i)]) && is.na(xparec[paste0("PMinB.", i)]) && is.na(xparec[paste0("Min.", i)]) && is.na(xparec[paste0("PMin.", i)])) {xparec[paste0("MinB.", i)] <- 0}
      if (is.na(xparec[paste0("MinE.", i)]) && is.na(xparec[paste0("PMinE.", i)]) && is.na(xparec[paste0("Min.", i)]) && is.na(xparec[paste0("PMin.", i)])) {xparec[paste0("MinE.", i)] <- 0}
      if (is.na(xparec[paste0("Flat.", i)])) {xparec[paste0("Flat.", i)] <- 0}
    }
  } else {
    
    if (is.na(xparec["MinB"]) && is.na(xparec["PMinB"]) && is.na(xparec["Min"]) && is.na(xparec["PMin"])) {xparec["MinB"] <- 0}
    if (is.na(xparec["MinE"]) && is.na(xparec["PMinE"]) && is.na(xparec["Min"]) && is.na(xparec["PMin"])) {xparec["MinE"] <- 0}
    xparec["Flat"] <- ifelse(is.na(xparec["Flat"]), 0, abs(xparec["Flat"]))
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
  
  if (index == 0) {
    if (!is.na(xparec["Length"])) {
      xparec["Begin"] <- xparec["Peak"] - xparec["Length"]
      xparec["End"] <- xparec["Peak"] + xparec["Length"]	
    }
    if (!is.na(xparec["LengthB"])) {
      xparec["Begin"] <- xparec["Peak"] - xparec["LengthB"]
    }
    if (!is.na(xparec["LengthE"])) {
      xparec["End"] <- xparec["Peak"] + xparec["LengthE"]	
    }
    
  } else {
    for (i in as.character(1:index)) {
      if ((!is.na(xparec["Length"])) & (is.na(xparec[paste0("Length.", i)]))) 
        xparec[paste0("Length.", i)] <- xparec["Length"]
      if ((!is.na(xparec["LengthB"])) & (is.na(xparec[paste0("LengthB.", i)]))) 
        xparec[paste0("LengthB.", i)] <- xparec["LengthB"]
      if ((!is.na(xparec["LengthE"])) & (is.na(xparec[paste0("LengthE.", i)]))) 
        xparec[paste0("LengthE.", i)] <- xparec["LengthE"]
      if ((!is.na(xparec["Peak"]))  & (is.na(xparec[paste0("Peak.", i)]))) 
        xparec[paste0("Peak.", i)] <- xparec["Peak"]
      
      if ((!is.na(xparec[paste0("Length.", i)])) & (is.na(xparec[paste0("Begin.", i)]))) {
        xparec[paste0("Begin.", i)] <- xparec[paste0("Peak.", i)] - xparec[paste0("Length.", i)]
      }
      if ((!is.na(xparec[paste0("Length.", i)])) & (is.na(xparec[paste0("End.", i)]))) {
        xparec[paste0("End.", i)] <- xparec[paste0("Peak.", i)] + xparec[paste0("Length.", i)]
      }
      if ((!is.na(xparec[paste0("LengthB.", i)])) & (is.na(xparec[paste0("Begin.", i)]))) {
        xparec[paste0("Begin.", i)] <- xparec[paste0("Peak.", i)] - xparec[paste0("LengthB.", i)]
      }
      if ((!is.na(xparec[paste0("LengthE.", i)])) & (is.na(xparec[paste0("End.", i)]))) {
        xparec[paste0("End.", i)] <- xparec[paste0("Peak.", i)]+xparec[paste0("LengthE.", i)]	
      }
    }
  }
  
  if (index == 0) {
    if (!is.na(xparec["PMinE"])) {xparec["MinE"]<-xparec["Max"]*xparec["PMinE"]/100}
    if (!is.na(xparec["PMinB"])) {xparec["MinB"]<-xparec["Max"]*xparec["PMinB"]/100}
    if (!is.na(xparec["PMin"])) {
      xparec["MinB"]<-xparec["Max"]*xparec["PMin"]/100
      xparec["MinE"]<-xparec["Max"]*xparec["PMin"]/100
    }
    if (!is.na(xparec["Min"])) {
      xparec["MinB"] <- xparec["Min"]
      xparec["MinE"] <- xparec["Min"]
    }
    
  } else {
    for (i in as.character(1:index)) {
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
    }
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
  
  if (index == 0) {
    xparec["PmoinsF"]<-xparec["Peak"]-(xparec["Flat"]/2)
    xparec["PplusF"]<-xparec["Peak"]+(xparec["Flat"]/2)
    
    xparec["PmoinsFB"]<-xparec["PmoinsF"]-xparec["Begin"]
    xparec["EPplusF"]<-xparec["End"]-xparec["PplusF"]
    
    xparec["MaxMinB"]<-xparec["Max"]-xparec["MinB"]
    xparec["MaxMinE"]<-xparec["Max"]-xparec["MinE"]
    
  } else {
    for (i in as.character(1:index)) {
      
      xparec[paste0("PmoinsF.", i)] <- xparec[paste0("Peak.", i)]-(xparec[paste0("Flat.", i)]/2)
      xparec[paste0("PplusF.", i)] <- xparec[paste0("Peak.", i)]+(xparec[paste0("Flat.", i)]/2)
      
      xparec[paste0("PmoinsFB.", i)] <- xparec[paste0("PmoinsF.", i)]-xparec[paste0("Begin.", i)]
      xparec[paste0("EPplusF.", i)] <- xparec[paste0("End.", i)]-xparec[paste0("PplusF.", i)]
      
      xparec[paste0("MaxMinB.", i)] <- xparec[paste0("Max.", i)]-xparec[paste0("MinB.", i)]
      xparec[paste0("MaxMinE.", i)] <- xparec[paste0("Max.", i)]-xparec[paste0("MinE.", i)]
      
    }
  }
  
  if (!is.null(season)) {
    xparec <- c(xparec[!grepl("\\.", names(xparec))], 
                xparec[grepl(paste0("\\.", season,"$"), names(xparec))])
    names(xparec) <- gsub("\\.[0-9]+$", "", names(xparec))
  }
  
  
  return(xparec)
  
}
