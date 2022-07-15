#' fitCF_MHmcmc_p generates set of parameters to be used with fitCF_MHmcmc()
#' @title Generate set of parameters to be used with fitCF_MHmcmc()
#' @author Marc Girondot
#' @return A matrix with the parameters
#' @param result An object obtained after a fitCF() fit
#' @param density Preset of density; can be dnorm or dunif
#' @param accept If TRUE, does not wait for user interaction
#' @family Model of Clutch Frequency
#' @description Interactive script used to generate set of parameters to be used with fitCF_MHmcmc().\cr
#' @examples 
#' \dontrun{
#' library("phenology")
#' data(MarineTurtles_2002)
#' ECFOCF_2002 <- TableECFOCF(MarineTurtles_2002)
#' 
#' # Paraetric model for clutch frequency
#' o_mu1p1_CFp <- fitCF(x = c(mu = 2.1653229641404539, 
#'                  sd = 1.1465246643327098, 
#'                  p = 0.25785366120357966), 
#'                  fixed.parameters=NULL, 
#'                  data=ECFOCF_2002, hessian = TRUE)
#'                            
#' pMCMC <- fitCF_MHmcmc_p(result=o_mu1p1_CFp, accept=TRUE)
#' fitCF_MCMC <- fitCF_MHmcmc(result = o_mu1p1_CFp, n.iter = 10000, 
#'                            parametersMCMC = pMCMC, n.chains = 1, n.adapt = 0, 
#'                            thin = 1, trace = FALSE)
#' plot(fitCF_MCMC, parameters="mu")
#' plot(fitCF_MCMC, parameters="sd")
#' plot(fitCF_MCMC, parameters="p", xlim=c(0, 0.5), breaks=seq(from=0, to=0.5, by=0.05))
#' plot(fitCF_MCMC, parameters="p", transform = invlogit, xlim=c(0, 1), 
#'      breaks=c(seq(from=0, to=1, by=0.05)))
#' }
#' @export

fitCF_MHmcmc_p <- function(result=stop("An output from fitCF() must be provided"), 
                           density="dunif", accept=FALSE) {
  
  if (!inherits(result, "ECFOCF")) {
    stop("An output from fitCF() must be provided")
  }
  
  
  # rownames(parametersMCMC)<-names(par)
  if (length(density) != length(result$par)) {
    density <- rep(density, length(result$par))[seq_along(result$par)]
  }
  
  parametersMCMC <- data.frame(Density=character(), Prior1=numeric(), Prior2=numeric(), 
                               SDProp=numeric(), Min=numeric(), Max=numeric(), Init=numeric(), stringsAsFactors = FALSE)
  
  for (indice.par in seq_along(result$par)) {
    par <- result$par[indice.par]
    SE <- result$SE[indice.par]
    if (is.na(SE)) SE <- par/2
    nm <- names(par)
    par <- unname(par)
    SE <- unname(SE)
    
    # J'ai une infinité de noms possibles !
    if (substr(nm, 1, 1)=="p") { 
      if (density[indice.par] == "dnorm") {
        parametersMCMC_ec <- data.frame(Density="dnorm", Prior1=par, Prior2=abs(SE), 
                                        SDProp=2, Min=-10*abs(par), 
                                        Max=10*abs(par), Init=par, 
                                        row.names = nm, stringsAsFactors = FALSE)
      } else {
        parametersMCMC_ec <- data.frame(Density="dunif", Prior1=-10*abs(par), Prior2=10*abs(par), 
                                        SDProp=2, Min=-10*abs(par), 
                                        Max=10*abs(par), Init=par, 
                                        row.names = nm, stringsAsFactors = FALSE)
      }
    } else {
      if (density[indice.par] == "dnorm") {
        parametersMCMC_ec <- data.frame(Density="dnorm", Prior1=par, Prior2=abs(SE), 
                                        SDProp=2, Min=0, 
                                        Max=10*par, Init=par, 
                                        row.names = nm, stringsAsFactors = FALSE)
      } else {
        parametersMCMC_ec <- data.frame(Density="dunif", Prior1=0, Prior2=10*par, 
                                        SDProp=2, Min=0, 
                                        Max=10*par, Init=par, 
                                        row.names = nm, stringsAsFactors = FALSE)
      }
    }
    parametersMCMC <- rbind(parametersMCMC, parametersMCMC_ec)
  }
  
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
            if (length(f)!=0) parameters[variable, "Prior1"] <- f
            cat(paste("Distribution of the prior, Maximum (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
            f<-scan(nmax=1, quiet=TRUE, what=character())
            if (length(f)!=0) parameters[variable, "Prior2"] <- f
            
          } else {
            
            if (density == "dnorm") {
              
              cat(paste("Distribution of the prior, Mean (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior1"] <- f
              cat(paste("Distribution of the prior, Standard deviation (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior2"] <- f
              
            } else {
              
              cat(paste("Distribution of the prior, value 1 (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior1"] <- f
              cat(paste("Distribution of the prior, value 2 (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior2"] <- f
              
            }
          }
          
          
          cat(paste("SD of new proposition (Enter for default ",parameters[variable, "SDProp"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "SDProp"] <- f
          cat(paste("Minimum for the parameter (default ",parameters[variable, "Min"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "Min"] <- f
          cat(paste("Maximum for the parameter (Enter for default ",parameters[variable, "Max"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "Max"] <- f
          cat(paste("Initial value (Enter for default ",parameters[variable, "Init"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "Init"] <- f
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
      warning(paste("Initial value for parameter ", rownames(parameters)[i], " was out of range; It is corrected. Check it.")) 
    }
  }
  
  
}
