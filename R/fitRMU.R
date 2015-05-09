#' fitRMU adjusts incomplete timeseries
#' @title Adjusts incomplete timeseries with various constraints.
#' @author Marc Girondot
#' @return Return a list with two components: df_final and result
#' @param data A data.frame with a column Year and two columns per rookery
#' @param method Can be Constant, Year-specific or Exponential
#' @param nbiter Maximum number of iterations before to show intermediate results
#' @param parameters  Parameters to fit
#' @param parametersfixed Parameters that are fixed
#' @param replicate.CI Number of replicates to estimate CI
#' @param colname.year Name of the column to be used as time index
#' @param RMU.name A dataframe with two columns indicating name of columns for mean and standard error for roockerys
#' @param optim Method of optim
#' @param control List for control parameters of optim
#' @description The data must be a data.frame with the first column being Year \cr
#' and two columns for each beach: the average and the se for the estimate.\cr
#' The correspondance between mean and se for each rookery are given in the RMU.name data.frame.\cr
#' @examples
#' \dontrun{
#' RMU.name.AtlanticW <- data.frame(mean=c("Yalimapo.French.Guiana", 
#'                                          "Galibi.Suriname", 
#'                                          "Irakumpapy.French.Guiana"), 
#'                                  se=c("se_Yalimapo.French.Guiana", 
#'                                       "se_Galibi.Suriname", 
#'                                       "se_Irakumpapy.French.Guiana"))
#' data.AtlanticW <- data.frame(Year=c(1990:2000), 
#'       Yalimapo.French.Guiana=c(2076, 2765, 2890, 2678, NA, 
#'                                6542, 5678, 1243, NA, 1566, 1566),
#'       se_Yalimapo.French.Guiana=c(123.2, 27.7, 62.5, 126, NA, 
#'                                  230, 129, 167, NA, 145, 20),
#'       Galibi.Suriname=c(276, 275, 290, NA, 267, 
#'                        542, 678, NA, 243, 156, 123),
#'       se_Galibi.Suriname=c(22.3, 34.2, 23.2, NA, 23.2, 
#'                            4.3, 2.3, NA, 10.3, 10.1, 8.9),
#'       Irakumpapy.French.Guiana=c(1076, 1765, 1390, 1678, NA, 
#'                                3542, 2678, 243, NA, 566, 566),
#'       se_Irakumpapy.French.Guiana=c(23.2, 29.7, 22.5, 226, NA, 
#'                                  130, 29, 67, NA, 15, 20))
#'                            
#' cst <- fitRMU(data=data.AtlanticW, RMU.name=RMU.name.AtlanticW, 
#'                colname.year="Year", method="Constant")
#' expo <- fitRMU(data=data.AtlanticW, RMU.name=RMU.name.AtlanticW, 
#'                colname.year="Year", method="Exponential")
#' YS <- fitRMU(data=data.AtlanticW, RMU.name=RMU.name.AtlanticW, 
#'                colname.year="Year", method="Year-specific")
#'                
#' compare_AIC(Constant=cst$result, Exponential=expo$result, 
#' YearSpecific=YS$result)
#' 
#' barplot_errbar(YS$mean, y.plus = YS$CI.0.95, 
#' y.minus = YS$CI.0.05, las=1, ylim=c(0, 0.7), 
#' main="Proportion of the different rookeries in the region")
#' 
#' plot_errbar(1990:2000, YS$result$par[paste0("T", 1990:2000)], 
#' bty="n", las=1, ylab="Nest number", xlab="Nesting season", 
#' errbar.y=2*YS$result$SE[paste0("T", 1990:2000)], ylim=c(0, 12000), 
#' main="Total number of nests in all the rookeries in the region")
#' }
#' @export


#################################################################################
# Ajustement
#################################################################################

fitRMU<-function(data=stop("data parameter must be provided"), 
                 method="Constant", 
                 nbiter=100, RMU.name=NULL, 
                 parameters=NULL, parametersfixed=NULL, optim="BFGS",
                 replicate.CI=1000, colname.year="Year", 
                 control=list(trace=1, REPORT=100, maxit=nbiter)) {
  
  # data=NA; method="Constant"; nbiter=100; RMU.name=NULL; parameters=NULL; parametersfixed=NULL; optim="BFGS"; replicate.CI=1000; colname.year="Year";control=list(trace=1, REPORT=100, maxit=nbiter)
  # data=data.AtlanticW; RMU.name=RMU.name.AtlanticW; colname.year="Year"; method="Constant"
  
  

  nm <- colnames(data)
  index.year <- which(nm==colname.year)
  index.mean <- match(RMU.name$mean, nm)
  index.se <- match(RMU.name$se, nm)
  
  index <- list(year=index.year, mean=index.mean, se=index.se)
  
  if (!all(c(index.year, index.mean, index.se) %in% seq_along(nm))) 
    stop("check the correspondance between names of columns and RMU.name or colname.year")
  
  # generate a vector of value, one per year
  # representing the total number of nests for all RMU for this year
  # The inital value is the mean value for the RMUs/years
  # ?????????
  # Nmean<-
  
  Tot <- NULL
  
  if (method=="Year-specific") {
    
    fun <- .fonctionfitYearspecific
    
    Tot <- rep(sum(colMeans(data[,index.mean], na.rm=TRUE), na.rm=TRUE), dim(data)[1])
    names(Tot) <- paste("T", data[,index.year], sep="")
  }
  
  if (method=="Constant") {
    # fait la moyenne colonne par colonne et fait la somme
    Tot <- sum(colMeans(data[,index.mean], na.rm=TRUE), na.rm=TRUE)
    names(Tot) <- "T"
    fun <- .fonctionfitTcst
  }
  
  
  if (method=="Exponential") {
    # fait la moyenne colonne par colonne et fait la somme
    Tot <- c(T=sum(colMeans(data[,index.mean], na.rm=TRUE), na.rm=TRUE), r=0)
    fun <- .fonctionfitTexp
  }
  

  if (!is.null(Tot)) {
    
    # generate a vector of values, one per RMU
    # representing the proportion of each RMU
    # The initial value is the mean value for each RMU
    
    p <- colMeans(data[,index.mean], na.rm=TRUE)
    last.index <- length(p)
    p<-(p/p[last.index])*100
    p<-p[-last.index]
    # attributes(p)$type<-rep("p", length(p))
    
    d<-data[,index.mean]
    
    
    x <- c(Tot, p, Size=mean(sd(d[!is.na(d)])))
    
    # ________________________________________________
    # Si un jeu de paramètre est indiqué, on l'utilise mais il n'est peut-être pas entier
    # il faut retirer ceux qui sont en commun
    # ________________________________________________
    if (any(!is.null(parameters))) {
      for(i in 1:length(parameters)) {
        x<-x[names(x)!=names(parameters[i])]
      }
      x<-c(x, parameters)
    }
    
    
    # ________________________________________________
    # S'il y a des paramètres fixes, on les retire du jeu de paramètres à ajuster
    # ________________________________________________
    
    if (!is.null(parametersfixed)) {
      for(i in 1:length(parametersfixed)) {
        x<-x[names(x)!=names(parametersfixed[i])]
      }
    }
    
    pfixed <- parametersfixed
    
    cpt<-0
    
    print("Initial values")
    print(x)
    
    repeat {
      result <- optim(x, fun, NULL, fixed=pfixed, data=data, index=index, method=optim, control=control, hessian=TRUE)
      
      if ((result$convergence==0)+(nbiter=0)) break
      cpt<-cpt+nbiter
      x<-result$par
      x[substr(names(x),1,1)=="p"]<-abs(x[substr(names(x),1,1)=="p"])
      x[substr(names(x),1,1)=="T"]<-abs(x[substr(names(x),1,1)=="T"])
      x["Size"]<-abs(x["Size"])
      print(paste("Iteration: ", cpt, " Convergence is not achieved -LnL=", format(result$value, digit=floor(log(result$value)/log(10))+3),sep=""))
      print(x)
    }
    x<-result$par
    print(paste("Convergence is achieved. -LnL=", format(result$value, digit=floor(log(result$value)/log(10))+3),sep=""))
    
    
    
    # Inverse the hessian matrix to get SE for each parameters
    mathessian=result$hessian;
    inversemathessian=solve(mathessian)
    res=sqrt(diag(inversemathessian))
    
    result$AIC <- 2*result$value+2*length(x)
    result$SE <- res
    
    print(paste("-Ln L=", result$value))
    print(paste("Parameters=", length(x)))
    print(paste("AIC=", result$AIC))
    
    
    print(x)
    
    pder <- 100
    names(pder) <- tail(names(data[,index.mean]), n=1L)
    seder <- 0
    names(seder) <- tail(names(data[,index.mean]), n=1L)
    
    p <- abs(c(result$par[names(result$par) %in% names(p)], pder))
    sep <- c(res[names(result$par) %in% names(p)], seder)
    
    randomp <- matrix(rep(p, replicate.CI), nrow=replicate.CI, byrow=TRUE)
    
    for(i in 1:length(p)) {
      randomp[,i]=abs(rnorm(replicate.CI, p[i], sep[i]))
    }
    
    for(j in 1:replicate.CI) {
      randomp[j,]=randomp[j,]/sum(randomp[j,])
    }
    
    sep2<-apply(randomp[,1:length(p)], 2, quantile, probs=c(0.05, 0.95))
    
    pmean<-p/sum(p)
    colnames(sep2) <- names(pmean)

    return(list(mean=pmean, CI.0.05=sep2[1,], CI.0.95=sep2[2,], result=result))
    
  } else {
    stop("Unknown method. Possible ones are: Year-specific, Constant, or Exponential")
  }
  
}


