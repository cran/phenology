#' map_phenology generates a likelihood map.
#' @title Generate a likelihood map varying Phi and Delta.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Display a likelihood map
#' @param data dataset generated with add_format
#' @param fixed.parameters Set of fixed parameters
#' @param fitted.parameters Set of parameters to be fitted
#' @param Phi Phi values to be analyzed
#' @param Delta Delta value to be analyzed
#' @param cofactors data.frame with a column Date and a column for each cofactor
#' @param add.cofactors Names of the column of parameter cofactors to use as a cofactor
#' @param zero If the theoretical nest number is under this value, this value wll be used
#' @param progressbar If FALSE, do not show the progress bar
#' @description This function generates a map of likelihood varying Phi and Delta.\cr
#' 	Parameters are the same than for the fit_phenology() function except for trace that is disabled.\cr
#' 	If Alpha, Beta or Tau are not indicated, Alpha and Tau are set to 0 and 1 and Beta is fitted.\cr
#' 	Only one set of Alpha, Beta, Tau, Phi and Delta are used for all timeseries present in data.\cr
#' 	Note that it is possible to fit or fixed Alpha[n], Beta[n], Tau[n], Phi[n] and Delta[n] with [n]=1 or 2 
#' 	and then it is possible to use this function to establish the likelihood map for a 
#' 	second or third sinusoids added to the global pattern.\cr
#' 	If Delta is not specified, it is estimated from Phi and the same precision as Phi is used.
#' @family Phenology model
#' @examples
#' library("phenology")
#' # Read a file with data
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot <- add_phenology(Gratiot, name = "Complete", 
#' 		reference = as.Date("2001-01-01"), format = "%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg <- par_init(data_Gratiot, fixed.parameters = NULL)
#' # Run the optimisation
#' \dontrun{
#' result_Gratiot <- fit_phenology(data = data_Gratiot, 
#' 		fitted.parameters = parg, fixed.parameters = NULL)
#' }
#' data(result_Gratiot)
#' # Extract the fitted parameters
#' parg1 <- extract_result(result_Gratiot)
#' # Add constant Alpha and Tau values 
#' # [day d amplitude=(Alpha+Nd*Beta)^Tau with Nd being the number of counts for day d]
#' pfixed <- c(parg1, Alpha=0, Tau=1)
#' pfixed <- pfixed[-which(names(pfixed)=="Theta")]
#' # The only fitted parameter will be Beta
#' parg2 <- c(Beta=0.5, parg1["Theta"])
#' # Generate a likelihood map 
#' # [default Phi=seq(from=0.1, to=20, length.out=100) but it is very long]
#' # Take care, it takes 20 hours ! The data map_Gratiot has the result
#' \dontrun{
#' library(phenology)
#' map_Gratiot <- map_phenology(data = data_Gratiot, 
#'                              Phi = seq(from=0.1, to=30, length.out=100), 
#' 		                          fitted.parameters = parg2, 
#' 		                          fixed.parameters = pfixed)
#' }
#' data(map_Gratiot)
#' # Plot the map
#' plot(map_Gratiot, col = heat.colors(128))
#' # Plot the min(-Ln L) for Phi varying at any delta value
#' plot_phi(map = map_Gratiot)
#' # Plot the min(-Ln L) for Delta varying with Phi equal to the value for maximum likelihood
#' plot_delta(map = map_Gratiot)
#' # Plot the min(-Ln L) for Delta varying with Phi the nearest to 15
#' plot_delta(map = map_Gratiot, Phi = 15)
#' @export

map_phenology <-
  function(data=NULL, fitted.parameters=NULL, fixed.parameters=NA, 
           Phi=seq(from=0.2,to=20, length.out=100), Delta=NULL, 
           progressbar=any(installed.packages()[, "Package"] == "pbapply"), 
           cofactors=NULL, add.cofactors=NULL, zero=1E-9) {
    
    # data=NULL; fitted.parameters=NULL; fixed.parameters=NA;Phi=seq(from=0.2,to=20, length.out=100); Delta=NULL; progressbar=TRUE; cofactors=NULL; add.cofactors=NULL; zero=1E-9
    mc.cores <- getOption("mc.cores", detectCores())
    forking <- getOption("forking", ifelse(.Platform$OS.type == "windows", FALSE, TRUE))
    
    if (is.null(fixed.parameters)) {fixed.parameters<-NA}
    if (is.null(fitted.parameters)) {fitted.parameters<-NA}
    
    # je varie Phi et Delta
    
    #create 2 vectors in form of numeric sequence, for Delta and Phi
    Phivalue <- Phi
    if (is.null(Delta)) {
      Deltavalue <- seq(from=0, to=max(Phivalue)/2, length.out=length(Phivalue)+1)
    } else {
      Deltavalue <- Delta
    }
    
    LPhi <- length(Phivalue)
    LDelta <- length(Deltavalue)
    
    # SET MATRIX
    # input <- matrix(data=NA, nrow=LPhi, ncol=LDelta)
    
    
    # si ni Alpha ni Beta ne sont à ajuster, je mets Beta
    if (is.na(fitted.parameters["Alpha"]) && is.na(fitted.parameters["Beta"])) {
      if (all(is.na(fitted.parameters))) {
        fitted.parameters <- c(Beta=0)
      } else {
        fitted.parameters <- c(fitted.parameters, Beta=0)
      }
    }
    
    # si Beta est à la fois fixe et à ajuster, je le retire des fixes
    if (!is.na(fitted.parameters["Beta"]) && !is.na(fixed.parameters["Beta"])) {fixed.parameters<-fixed.parameters[!names(fixed.parameters)=="Beta"]}
    
    # je vérifie que Alpha, Beta et Tau apparaissent bien au moins une fois, sinon je les mets en fixe
    xpar <- c(fitted.parameters, fixed.parameters)
    if (is.na(xpar["Alpha"])) {fixed.parameters<-c(fixed.parameters, Alpha=0)}
    if (is.na(xpar["Beta"])) {fixed.parameters<-c(fixed.parameters, Beta=0)}
    if (is.na(xpar["Tau"])) {fixed.parameters<-c(fixed.parameters, Tau=1)}
    
    # Si Phi ou Delta sont indiqués en paramètres à ajuster, je les retire
    fitted.parameters <- fitted.parameters[!names(fitted.parameters)=="Phi"]
    fitted.parameters <- fitted.parameters[!names(fitted.parameters)=="Delta"]
    
    # mais je les mets en paramètres fixes
    if (is.na(fixed.parameters["Phi"])) {fixed.parameters <- c(fixed.parameters, Phi=0)}
    if (is.na(fixed.parameters["Delta"])) {fixed.parameters <- c(fixed.parameters, Delta=0)}
    
    
    # if (progressbar) pb <- txtProgressBar(min=0, max=LDelta, style=3)
    
    #parpre1<-fitted.parameters
    
    if ((progressbar) & (any(installed.packages()[, "Package"] == "pbapply"))) {
      CEG <- list(expr=expression(library("phenology"), 
                                  requireNamespace("pbapply", quietly=TRUE)))
      # requireNamespace("pbapply", quietly=TRUE)
    } else {
      CEG <- list(expr=expression(library("phenology")))
    }
    
    pt <- list(data=data, fixed=fixed.parameters, out=TRUE, 
               cofactors=cofactors, 
               add.cofactors=add.cofactors, zero=zero, 
               store.intermediate=FALSE, 
               file.intermediate="")
    
    # ptx<<-pt
    
    outma <- universalmclapply(1:(LDelta*LPhi), FUN=function(ma) {
      
      #FILLING MATRIX
      # for(j in 1:LDelta) {
      

        
        map <- ma - 1
        
        j <- (map %/% LPhi) + 1 # lignes
        i <- (map %% LPhi) + 1 # colonnes
        
        XDelta <- Deltavalue[j]
        #   for(i in 1:LPhi) {
        
        XPhi <- Phivalue[i]
        
        # pt$fixed["Delta"] <- XDelta
        # return(XDelta)
        
        if (XDelta >= XPhi/2) {
          out <- NA
        } else {
          
          pt$fixed["Delta"] <- XDelta
          pt$fixed["Phi"] <- XPhi
          
          #	assign("fixed", fixed.parameters, envir=as.environment(.phenology.env))
          
          par <- fitted.parameters
          
          Lpre <- 0
          
          
          repeat {
            resul <- optim(par, getFromNamespace(".Lnegbin", ns="phenology"), 
                           pt=pt, 
                           method="BFGS", 
                           control=list(trace=0, REPORT=1, maxit=5000), 
                           hessian=FALSE)
            if ((resul$convergence == 0) & (abs(resul$value-Lpre) < 1E-5)) break
            par <- resul$par
            Lpre <- resul$value
          }
          
          out <- resul$value
          
        }
        return(out)
        
    }, clusterExport = list(pt=pt, 
                            fixed.parameters = fixed.parameters,
                            LPhi = LPhi, 
                            Deltavalue = Deltavalue, 
                            Phivalue = Phivalue), 
    clusterEvalQ = CEG, 
    progressbar=progressbar, forking=forking, mc.cores=mc.cores)
    
    input <- matrix(unlist(outma), nrow=LPhi, ncol=LDelta, byrow=FALSE)
    rownames(input) <- Phivalue
    colnames(input) <- Deltavalue
    
    Dv <- as.vector(input)
    pos <- which(x = (input == min(input, na.rm = TRUE)), arr.ind = TRUE)
    i0 <- pos[1, 1]
    j0 <- pos[1, 2]
    print(paste("The minimum -Ln likelihood is ", input[i0, j0], sep=""))
    print(paste("For Phi=",Phivalue[i0],sep=""))
    print(paste("And Delta=",Deltavalue[j0],sep=""))
    
    outputmap <- list(input=input, Phi=Phivalue, Delta=Deltavalue, 
                      Fitted.parameters=names(fitted.parameters), 
                      Fixed.parameters=fixed.parameters, 
                      Data=names(data))
    
    outputmap <- addS3Class(outputmap, "phenologymap")
    return(outputmap)
    
  }
