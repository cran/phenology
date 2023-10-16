#' phenology_MHmcmc_p generates set of parameters to be used with phenology_MHmcmc()
#' @title Generates set of parameters to be used with phenology_MHmcmc()
#' @author Marc Girondot
#' @return A matrix with the parameters
#' @param result An object obtained after a fit_phenology() fit
#' @param default.density The default density, "dnorm" or "dunif'
#' @param accept If TRUE, does not wait for use interaction
#' @description Interactive script used to generate set of parameters to be used with phenology_MHmcmc().
#' @family Phenology model
#' @examples 
#' \dontrun{
#' library(phenology)
#' data(Gratiot)
#' # Generate a formatted list named data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#'   	reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, fixed.parameters=NULL)
#' # Run the optimisation
#' result_Gratiot<-fit_phenology(data=data_Gratiot, 
#' 		fitted.parameters=parg, fixed.parameters=NULL)
#' # Generate set of priors for Bayesian analysis
#' pmcmc <- phenology_MHmcmc_p(result_Gratiot, accept = TRUE)
#' result_Gratiot_mcmc <- phenology_MHmcmc(result = result_Gratiot, n.iter = 10000, 
#' parametersMCMC = pmcmc, n.chains = 1, n.adapt = 0, thin = 1, trace = FALSE)
#' # Get standard error of parameters
#' summary(result_Gratiot_mcmc)
#' # Make diagnostics of the mcmc results using coda package
#' mcmc <- as.mcmc(result_Gratiot_mcmc)
#' require(coda)
#' heidel.diag(mcmc)
#' raftery.diag(mcmc)
#' autocorr.diag(mcmc)
#' acf(mcmc[[1]][,"LengthB"], lag.max=200, bty="n", las=1)
#' acf(mcmc[[1]][,"Max_Gratiot"], lag.max=50, bty="n", las=1)
#' batchSE(mcmc, batchSize=100)
#' # The batch standard error procedure is usually thought to 
#' # be not as accurate as the time series methods used in summary
#' summary(mcmc)$statistics[,"Time-series SE"]
#' plot(result_Gratiot_mcmc, parameters=3, las=1, xlim=c(-10, 300))
#' }
#' @export

phenology_MHmcmc_p<-function(result=stop("An output from fit_phenology() must be provided"), 
                             default.density="dunif", 
                             accept=FALSE) {
  
  if (!inherits(result, "phenology")) {
    stop("An output of fit_phenology() must be provided")
  }
  
  # d'abord je sors les paramètres à utiliser
  
  par <- result$par
  
  # "Peak"
  pe <- ifelse(is.na(par["Peak"]), 180, par["Peak"])
  if (default.density == "dunif") {
    Peak <-  c("dunif", 0, max(c(365, pe+100)), 5, 0, max(c(365, pe+100)), pe)
  } else {
    Peak <-  c("dnorm", pe, pe/2, 5, 0, max(c(365, pe+100)), pe)
  }
  
  # "Flat"
  pe <- ifelse(is.na(par["Flat"]), 5, par["Flat"])
  if (default.density == "dunif") {
    Flat <- c("dunif", 0, max(c(50, pe+20)), 2, 0, max(c(50, pe+20)), pe)
  } else {
    Flat <- c("dnorm", pe, pe/2, 2, 0, max(c(50, pe+20)), pe)
  }
  
  # "Begin"
  pe <- ifelse(is.na(par["Begin"]), 100, par["Begin"])
  if (default.density == "dunif") {
    Begin <- c("dunif", 0, max(c(365, pe+50)), 20, 0, max(c(365, pe+50)), pe)
  } else {
    Begin <- c("dnorm", pe, pe/2, 20, 0, max(c(365, pe+50)), pe)
  }
  
  # "End"
  pe <- ifelse(is.na(par["End"]), 220, par["End"])
  if (default.density == "dunif") {
    End <- c("dunif", 0, max(c(365, pe+50)), 20, 0, max(c(365, pe+50)), pe)
  } else {
    End <- c("dnorm", pe, pe/2, 20, 0, max(c(365, pe+50)), pe)
  }
  
  # "Length"
  pe <- ifelse(is.na(par["Length"]), 100, par["Length"])
  if (default.density == "dunif") {
    Length <- c("dunif", 0, max(c(200, pe+50)), 20, 0, max(c(200, pe+50)), pe)
  } else {
    Length <- c("dnorm", pe, pe/2, 20, 0, max(c(200, pe+50)), pe)
  }
  
  # "LengthE"
  pe <- ifelse(is.na(par["LengthE"]), 100, par["LengthE"])
  if (default.density == "dunif") {
    LengthE <- c("dunif", 0, max(c(200, pe+50)), 20, 0, max(c(200, pe+50)), pe)
  } else {
    LengthE <- c("dnorm", pe, pe/2, 20, 0, max(c(200, pe+50)), pe)
  }
  
  # "LengthB"
  pe <- ifelse(is.na(par["LengthB"]), 100, par["LengthB"])
  if (default.density == "dunif") {
    LengthB <- c("dunif", 0, max(c(200, pe+50)), 20, 0, max(c(200, pe+50)), pe)
  } else {
    LengthB <- c("dnorm", pe, pe/2, 20, 0, max(c(200, pe+50)), pe)
  }
  
  # "Max"
  pe <- ifelse(is.na(par["Max"]), 50, par["Max"])
  if (default.density == "dunif") {
    Max <-c("dunif", 0, max(c(200, pe+50)), 0.4, 0, max(c(200, pe+50)), pe)
  } else {
    Max <-c("dnorm", pe, pe/2, 0.4, 0, max(c(200, pe+50)), pe)
  }
  
  # "PMin"
  pe <- ifelse(is.na(par["PMin"]), 100, par["PMin"])
  if (default.density == "dunif") {
    PMin <- c("dunif", 0, max(c(10, pe+10)), 2, 0, max(c(10, pe+10)), pe)
  } else {
    PMin <- c("dnorm", pe, pe/2, 0.5, 0, max(c(10, pe+10)), pe)
  }
  
  # "Min"
  pe <- ifelse(is.na(par["Min"]), 100, par["Min"])
  if (default.density == "dunif") {
    Min <- c("dunif", 0, max(c(5, pe+5)), 0.5, 0, max(c(5, pe+5)), pe)
  } else {
    Min <- c("dnorm", pe, pe/2, 0.1, 0, max(c(5, pe+5)), pe)
  }
  
  # "PMinE"
  pe <- ifelse(is.na(par["PMinE"]), 100, par["PMinE"])
  if (default.density == "dunif") {
    PMinE <- c("dunif", 0, max(c(10, pe+10)), 2, 0, max(c(10, pe+10)), pe)
  } else {
    PMinE <- c("dnorm", pe, pe/2, 0.5, 0, max(c(10, pe+10)), pe)
  }
  
  # "MinE"
  pe <- ifelse(is.na(par["MinE"]), 100, par["MinE"])
  if (default.density == "dunif") {
    MinE <- c("dunif", 0, max(c(5, pe+5)), 0.5, 0, max(c(5, pe+5)), pe)
  } else {
    MinE <- c("dnorm", pe, pe/2, 0.1, 0, max(c(5, pe+5)), pe)
  }
  
  # "PMinB"
  pe <- ifelse(is.na(par["PMinB"]), 100, par["PMinB"])
  if (default.density == "dunif") {
    PMinB <- c("dunif", 0, max(c(10, pe+10)), 2, 0, max(c(10, pe+10)), pe)
  } else {
    PMinB <- c("dnorm", pe, pe/2, 0.5, 0, max(c(10, pe+10)), pe)
  }
  
  # "MinB"
  pe <- ifelse(is.na(par["MinB"]), 100, par["MinB"])
  if (default.density == "dunif") {
    MinB <- c("dunif", 0, max(c(5, pe+5)), 0.5, 0, max(c(5, pe+5)), pe)
  } else {
    MinB <- c("dnorm", pe, pe/2, 0.1, 0, max(c(5, pe+5)), pe)
  }
  
  # "Phi"
  pe <- ifelse(is.na(par["Phi"]), 30, par["Phi"])
  if (default.density == "dunif") {
    Phi <- c("dunif", 0, max(c(50, pe+50)), 3, 0, max(c(50, pe+50)), pe)
  } else {
    Phi <- c("dnorm", pe, pe/2, 3, 0, max(c(50, pe+50)), pe)
  }
  
  # "Delta"
  pe <- ifelse(is.na(par["Delta"]), 30, par["Delta"])
  if (default.density == "dunif") {
    Delta <- c("dunif", 0, max(c(50, pe+50)), 3, 0, max(c(50, pe+50)), pe)
  } else {
    Delta <- c("dnorm", pe, pe/2, 3, 0, max(c(50, pe+50)), pe)
  }
  
  # "Alpha"
  pe <- ifelse(is.na(par["Alpha"]), 30, par["Alpha"])
  if (default.density == "dunif") {
    Alpha <- c("dunif", 0, max(c(50, pe+50)), 3, 0, max(c(50, pe+50)), pe)
  } else {
    Alpha <- c("dnorm", pe, pe/2, 3, 0, max(c(50, pe+50)), pe)
  }
  
  # "Beta"
  pe <- ifelse(is.na(par["Beta"]), 30, par["Beta"])
  if (default.density == "dunif") {
    Beta <- c("dunif", 0, max(c(50, pe+50)), 3, 0, max(c(50, pe+50)), pe)
  } else {
    Beta <- c("dnorm", pe, pe/2, 3, 0, max(c(50, pe+50)), pe)
  }
  
  # "Tau"
  pe <- ifelse(is.na(par["Tau"]), 2, par["Tau"])
  if (default.density == "dunif") {
    Tau <- c("dunif", 0, max(c(5, pe+5)), 0.5, 0, max(c(5, pe+5)), pe)
  } else {
    Tau <- c("dnorm", pe, pe/2, 0.5, 0, max(c(5, pe+5)), pe)
  }
  
  # "Phi1"
  pe <- ifelse(is.na(par["Phi1"]), 30, par["Phi1"])
  if (default.density == "dunif") {
    Phi1 <- c("dunif", 0, max(c(50, pe+50)), 3, 0, max(c(50, pe+50)), pe)
  } else {
    Phi1 <- c("dnorm", pe, pe/2, 3, 0, max(c(50, pe+50)), pe)
  }
  
  # "Delta1"
  pe <- ifelse(is.na(par["Delta1"]), 30, par["Delta1"])
  if (default.density == "dunif") {
    Delta1 <- c("dunif", 0, max(c(50, pe+50)), 3, 0, max(c(50, pe+50)), pe)
  } else {
    Delta1 <- c("dnorm", pe, pe/2, 3, 0, max(c(50, pe+50)), pe)
  }
  
  # "Alpha1"
  pe <- ifelse(is.na(par["Alpha1"]), 30, par["Alpha1"])
  if (default.density == "dunif") {
    Alpha1 <- c("dunif", 0, max(c(50, pe+50)), 3, 0, max(c(50, pe+50)), pe)
  } else {
    Alpha1 <- c("dnorm", pe, pe/2, 3, 0, max(c(50, pe+50)), pe)
  }
  
  # "Beta1"
  pe <- ifelse(is.na(par["Beta1"]), 30, par["Beta1"])
  if (default.density == "dunif") {
    Beta1 <- c("dunif", 0, max(c(50, pe+50)), 3, 0, max(c(50, pe+50)), pe)
  } else {
    Beta1 <- c("dnorm", pe, pe/2, 3, 0, max(c(50, pe+50)), pe)
  }
  
  # "Tau1"
  pe <- ifelse(is.na(par["Tau1"]), 2, par["Tau1"])
  if (default.density == "dunif") {
    Tau1 <- c("dunif", 0, max(c(5, pe+5)), 0.5, 0, max(c(5, pe+5)), pe)
  } else {
    Tau1 <- c("dnorm", pe, pe/2, 0.5, 0, max(c(5, pe+5)), pe)
  }
  
  # "Phi2"
  pe <- ifelse(is.na(par["Phi2"]), 30, par["Phi2"])
  if (default.density == "dunif") {
    Phi2 <- c("dunif", 0, max(c(50, pe+50)), 3, 0, max(c(50, pe+50)), pe)
  } else {
    Phi2 <- c("dnorm", pe, pe/2, 3, 0, max(c(50, pe+50)), pe)
  }
  
  # "Delta2"
  pe <- ifelse(is.na(par["Delta2"]), 30, par["Delta2"])
  if (default.density == "dunif") {
    Delta2 <- c("dunif", 0, max(c(50, pe+50)), 3, 0, max(c(50, pe+50)), pe)
  } else {
    Delta2 <- c("dnorm", pe, pe/2, 3, 0, max(c(50, pe+50)), pe)
  }
  
  # "Alpha2"
  pe <- ifelse(is.na(par["Alpha2"]), 30, par["Alpha2"])
  if (default.density == "dunif") {
    Alpha2 <- c("dunif", 0, max(c(50, pe+50)), 3, 0, max(c(50, pe+50)), pe)
  } else {
    Alpha2 <- c("dnorm", pe, pe/2, 3, 0, max(c(50, pe+50)), pe)
  }
  
  # "Beta2"
  pe <- ifelse(is.na(par["Beta2"]), 30, par["Beta2"])
  if (default.density == "dunif") {
    Beta2 <- c("dunif", 0, max(c(50, pe+50)), 3, 0, max(c(50, pe+50)), pe)
  } else {
    Beta2 <- c("dnorm", pe, pe/2, 3, 0, max(c(50, pe+50)), pe)
  }
  
  # "Tau2"
  pe <- ifelse(is.na(par["Tau2"]), 2, par["Tau2"])
  if (default.density == "dunif") {
    Tau2 <- c("dunif", 0, max(c(5, pe+5)), 0.5, 0, max(c(5, pe+5)), pe)
  } else {
    Tau2 <- c("dnorm", pe, pe/2, 0.5, 0, max(c(5, pe+5)), pe)
  }
  
  # "Theta"
  pe <- ifelse(is.na(par["Theta"]), 5, par["Theta"])
  if (default.density == "dunif") {
    Theta <- c("dunif", 1E-6, max(c(10, pe+5)), 0.2, 1E-6, max(c(10, pe+5)), pe)
  } else {
    Theta <- c("dnorm", pe, pe/2, 0.2, 1E-6, max(c(10, pe+5)), pe)
  }
  
  # "alpha"
  pe <- ifelse(is.na(par["alpha"]), 5, par["alpha"])
  if (default.density == "dunif") {
    alpha <-c("dunif", 0, max(c(200, pe+50)), 0.4, 0, max(c(200, pe+50)), pe)
  } else {
    alpha <-c("dnorm", pe, pe/2, 0.4, 0, max(c(200, pe+50)), pe)
  }
  
  # "tp"
  pe <- ifelse(is.na(par["tp"]), 180, par["tp"])
  if (default.density == "dunif") {
    tp <-  c("dunif", 0, max(c(365, pe+100)), 5, 0, max(c(365, pe+100)), pe)
  } else {
    tp <-  c("dnorm", pe, pe/2, 5, 0, max(c(365, pe+100)), pe)
  }
  
  # "tf"
  pe <- ifelse(is.na(par["tf"]), 40, par["tf"])
  if (default.density == "dunif") {
    tf <-  c("dunif", 0, max(c(40, pe+50)), 20, 0, max(c(40, pe+50)), pe)
  } else {
    tf <-  c("dnorm", pe, pe/2, 20, 0, max(c(40, pe+50)), pe)
  }
  
  # "s1"
  pe <- ifelse(is.na(par["s1"]), 100, par["s1"])
  if (default.density == "dunif") {
    s1 <- c("dunif", 0, max(c(200, pe+50)), 20, 0, max(c(200, pe+50)), pe)
  } else {
    s1 <- c("dnorm", pe, pe/2, 20, 0, max(c(200, pe+50)), pe)
  }
  
  # "s2"
  pe <- ifelse(is.na(par["s2"]), 100, par["s2"])
  if (default.density == "dunif") {
    s2 <- c("dunif", 0, max(c(200, pe+50)), 20, 0, max(c(200, pe+50)), pe)
  } else {
    s2 <- c("dnorm", pe, pe/2, 20, 0, max(c(200, pe+50)), pe)
  }
  
  # "sr"
  pe <- ifelse(is.na(par["sr"]), 100, par["sr"])
  if (default.density == "dunif") {
    sr <- c("dunif", 0, max(c(200, pe+50)), 20, 0, max(c(200, pe+50)), pe)
  } else {
    sr <- c("dnorm", pe, pe/2, 20, 0, max(c(200, pe+50)), pe)
  }
  
  
  priors <- list(Peak, Flat, Begin, End, Length, LengthE, LengthB, 
                 Length, Max, PMin, Min, PMinE, MinE, PMinB, MinB, Phi, Delta, Alpha, 
                 Beta, Tau, Phi1, Delta1, Alpha1, Beta1, Tau1, Phi2, Delta2, 
                 Alpha2, Beta2, Tau2, Theta, alpha, tp, tf, s1, s2, sr)
  
  names(priors) <- c("Peak", "Flat", "Begin", "End", "Length", 
                     "LengthE", "LengthB", "Length", "Max", "PMin", "Min", "PMinE", "MinE", 
                     "PMinB", "MinB", "Phi", "Delta", "Alpha", "Beta", "Tau", "Phi1", 
                     "Delta1", "Alpha1", "Beta1", "Tau1", "Phi2", "Delta2", "Alpha2", 
                     "Beta2", "Tau2", "Theta", "alpha", "tp", "tf", "s1", "s2", "sr")
  
  for (i in seq_along(par)) {
    
    if (substr(names(par[i]), 1, 4)=="Max_") {
      pe <- ifelse(is.na(par[i]), 50, par[i])
      if (default.density == "dunif") {
        priors <- c(priors, list(c("dunif", 0, max(c(200, pe+50)), 2, 0, max(c(200, pe+50)), pe)))
      } else {
        if (pe != 0) {
          priors <- c(priors, list(c("dnorm", pe, pe/2, pe/10, 0, max(c(200, pe+50)), pe)))
        } else {
          priors <- c(priors, list(c("dnorm", pe, 10, 2, 0, max(c(200, pe+50)), pe))) 
        }
      }
      names(priors)[length(priors)] <- names(par[i])
    }
    
    if (substr(names(par[i]), 1, 4)=="Min_") {
      pe <- ifelse(is.na(par[i]), 5, par[i])
      if (default.density == "dunif") {
        priors <- c(priors, list(c("dunif", 0, max(c(5, pe+5)), 2, 0, max(c(5, pe+5)), pe)))
      } else {
        priors <- c(priors, list(c("dnorm", pe, pe/2, pe/10, 0, max(c(5, pe+5)), pe)))
      }
      names(priors)[length(priors)] <- names(par[i])
    }
    
    if (substr(names(par[i]), 1, 5)=="MinE_") {
      pe <- ifelse(is.na(par[i]), 5, par[i])
      if (default.density == "dunif") {
        priors <- c(priors, list(c("dunif", 0, max(c(5, pe+5)), 2, 0, max(c(5, pe+5)), pe)))
      } else {
        priors <- c(priors, list(c("dnorm", pe, pe/2, pe/10, 0, max(c(5, pe+5)), pe)))
      }
      names(priors)[length(priors)] <- names(par[i])
    }
    
    if (substr(names(par[i]), 1, 5)=="MinB_") {
      pe <- ifelse(is.na(par[i]), 5, par[i])
      if (default.density == "dunif") {
        priors <- c(priors, list(c("dunif", 0, max(c(5, pe+5)), 2, 0, max(c(5, pe+5)), pe)))
      } else {
        priors <- c(priors, list(c("dnorm", pe, pe/2, pe/10, 0, max(c(5, pe+5)), pe)))
      }
      names(priors)[length(priors)] <- names(par[i])
    }
    if (substr(names(par[i]), 1, 5)=="Peak_") {
      pe <- ifelse(is.na(par[i]), 180, par[i])
      if (default.density == "dunif") {
        priors <-  c(priors, list(c("dunif", 0, max(c(365, pe+100)), 5, 0, max(c(365, pe+100)), pe)))
      } else {
        priors <-  c(priors, list(c("dnorm", pe, pe/2, 5, 0, max(c(365, pe+100)), pe)))
      }
      names(priors)[length(priors)] <- names(par[i])
    }
    if (substr(names(par[i]), 1, 7)=="Length_") {
      pe <- ifelse(is.na(par[i]), 100, par[i])
      if (default.density == "dunif") {
        priors <- c(priors, list(c("dunif", 0, max(c(200, pe+50)), 20, 0, max(c(200, pe+50)), pe)))
      } else {
        priors <- c(priors, list(c("dnorm", pe, pe/2, 20, 0, max(c(200, pe+50)), pe)))
      }
      names(priors)[length(priors)] <- names(par[i])
    }
    if (substr(names(par[i]), 1, 8)=="LengthB_") {
      pe <- ifelse(is.na(par[i]), 100, par[i])
      if (default.density == "dunif") {
        priors <- c(priors, list(c("dunif", 0, max(c(200, pe+50)), 20, 0, max(c(200, pe+50)), pe)))
      } else {
        priors <- c(priors, list(c("dnorm", pe, pe/2, 20, 0, max(c(200, pe+50)), pe)))
      }
      names(priors)[length(priors)] <- names(par[i])
    }
    if (substr(names(par[i]), 1, 8)=="LengthE_") {
      pe <- ifelse(is.na(par[i]), 100, par[i])
      if (default.density == "dunif") {
        priors <- c(priors, list(c("dunif", 0, max(c(200, pe+50)), 20, 0, max(c(200, pe+50)), pe)))
      } else {
        priors <- c(priors, list(c("dnorm", pe, pe/2, 20, 0, max(c(200, pe+50)), pe)))
      }
      names(priors)[length(priors)] <- names(par[i])
    }
    if (substr(names(par[i]), 1, 6)=="Begin_") {
      pe <- ifelse(is.na(par[i]), 100, par[i])
      if (default.density == "dunif") {
        priors <- c(priors, list(c("dunif", 0, max(c(365, pe+50)), 20, 0, max(c(365, pe+50)), pe)))
      } else {
        priors <- c(priors, list(c("dnorm", pe, pe/2, 20, 0, max(c(365, pe+50)), pe)))
      }
      names(priors)[length(priors)] <- names(par[i])
    }
    if (substr(names(par[i]), 1, 4)=="End_") {
      pe <- ifelse(is.na(par[i]), 100, par[i])
      if (default.density == "dunif") {
        priors <- c(priors, list(c("dunif", 0, max(c(365, pe+50)), 20, 0, max(c(365, pe+50)), pe)))
      } else {
        priors <- c(priors, list(c("dnorm", pe, pe/2, 20, 0, max(c(365, pe+50)), pe)))
      }
      names(priors)[length(priors)] <- names(par[i])
    }
    if (substr(names(par[i]), 1, 6)=="Theta_") {
      pe <- ifelse(is.na(par[i]), 1, par[i])
      if (default.density == "dunif") {
        priors <- c(priors, list(c("dunif", 1E-6, max(c(10, pe+5)), 2, 1E-6, max(c(10, pe+5)), pe)))
      } else {
        priors <- c(priors, list(c("dnorm", pe, pe/2, 2, 1E-6, max(c(10, pe+5)), pe)))
      }
      names(priors)[length(priors)] <- names(par[i])
    }
    
  }
  
  
  prencours <- NULL
  perror <- NULL
  for (i in seq_along(par)) {
    if (!is.null(priors[[names(par)[i]]])) {
    prencours <- c(prencours, priors[[names(par)[i]]])
    } else {
      perror <- c(perror, names(par)[i])
    }
  }
  
  if (length(perror) != 0) {
    stop(paste0(c("The parameters", perror, "are not managed by this function; you must build the priors by yourself."), collapse = ", "))
  }
  
  parametersMCMC <- matrix(prencours, ncol=7, byrow=T)
  colnames(parametersMCMC) <- c("Density", "Prior1", "Prior2", "SDProp", "Min", "Max", "Init")
  rownames(parametersMCMC)<- names(par)
  parametersMCMC <- as.data.frame(parametersMCMC, stringsAsFactors = FALSE)
  
  for (i in 2:7)
    parametersMCMC[,i] <- as.numeric(parametersMCMC[,i])
  
  
  parameters <- parametersMCMC
  
  if (accept) {
    return(parameters)
  } else {
    
    repeat {
      
      cat("Proposition:\n")
      print(parameters)
      cat("Name of the parameter to change or Enter to quit:\n")
      f<-scan(nmax=1, quiet=TRUE, what=character())
      
      if (length(f)==0) f <- "q"
      
      if (f=="q") {
        return(parameters)
        
      } else {
        
        variable <- which(f==names(par))
        if (length(variable)==0) {
          cat("The parameter does not exist:\n")
        } else {
          print(variable)
          cat(paste("Change for the parameter ",names(par)[variable],":\n",sep=""))
          
          cat(paste("Distribution of the prior (Enter for default ",parameters[variable, "Density"], "):", sep=""))
          density<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(density)!=0) { parameters[variable, "Density"] <- density } else { density <- parameters[variable, "Density"] }
          
          if (density == "dunif") {
            
            cat(paste("Distribution of the prior, Minimum (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
            f<-scan(nmax=1, quiet=TRUE, what=character())
            if (length(f)!=0) parameters[variable, "Prior1"] <- as.numeric(f)
            cat(paste("Distribution of the prior, Maximum (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
            f<-scan(nmax=1, quiet=TRUE, what=character())
            if (length(f)!=0) parameters[variable, "Prior2"] <- as.numeric(f)
            
          } else {
            
            if (density == "dnorm") {
              
              cat(paste("Distribution of the prior, Mean (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior1"] <- as.numeric(f)
              cat(paste("Distribution of the prior, Standard deviation (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior2"] <- as.numeric(f)
              
            } else {
              
              cat(paste("Distribution of the prior, value 1 (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior1"] <- as.numeric(f)
              cat(paste("Distribution of the prior, value 2 (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior2"] <- as.numeric(f)
              
            }
          }
          
          
          cat(paste("SD of new proposition (Enter for default ",parameters[variable, "SDProp"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "SDProp"] <- as.numeric(f)
          cat(paste("Minimum for the parameter (default ",parameters[variable, "Min"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "Min"] <- as.numeric(f)
          cat(paste("Maximum for the parameter (Enter for default ",parameters[variable, "Max"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "Max"] <- as.numeric(f)
          cat(paste("Initial value (Enter for default ",parameters[variable, "Init"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "Init"] <- as.numeric(f)
        }
        
      }
      
    }
    
  }
  
  for (i in 1:nrow(parameters)) {
    if (parameters[i, "Density"]=="dunif") {
      mn <- max(as.numeric(parameters[i, "Prior1"]), as.numeric(parameters[i, "Min"]))    
    } else {
      mn <- as.numeric(parameters[i, "Min"])
      mx <- as.numeric(parameters[i, "Max"])
    }  
    if (findInterval(as.numeric(parameters[i, "Init"]), c(mn, mx)) != 1) {
      parameters[i, "Init"] <- as.character(mn+(mx-mn)/2)
      warning(paste0("Initial value for parameter ", rownames(parameters)[i], " was out of range; It has been corrected. Check it.")) 
    }
  }
  
  
}
