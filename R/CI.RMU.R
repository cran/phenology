#' CI.RMU calculates the confidence interval of the results of fitRMU()
#' @title Calculate the confidence interval of the results of fitRMU()
#' @author Marc Girondot
#' @return Return a list with Total, Proportions, and Numbers
#' @param result A result of fitRMu()
#' @param resultMCMC A resuts of fitRMU_MHmcmc()
#' @param chain Number of MCMC chain to be used
#' @param regularThin If TRUE, use regular thin for MCMC
#' @param replicate.CI Number of replicates
#' @param silent If TRUE does not display anything
#' @family Fill gaps in RMU
#' @description The data must be a data.frame with the first column being years \cr
#' and two columns for each beach: the average and the se for the estimate.\cr
#' The correspondence between mean, se and density for each rookery are given in the RMU.names data.frame.\cr
#' This data.frame must have a column named mean, another named se and a third named density. If 
#' no sd column exists, no sd will be considered for the series and is no density column exists, it 
#' will be considered as being "dnorm".\cr
#' In the result list, the mean proportions for each rookeries are in $proportions, $proportions.CI.0.05 and $proportions.CI.0.95.\cr
#' The names of beach columns must not begin by T_, SD_, a0_, a1_ or a2_ and cannot be r.\cr
#' A RMU is the acronyme for Regional Managment Unit. See:\cr
#' Wallace, B.P., DiMatteo, A.D., Hurley, B.J., Finkbeiner, E.M., Bolten, A.B., 
#' Chaloupka, M.Y., Hutchinson, B.J., Abreu-Grobois, F.A., Amorocho, D., Bjorndal, K.A., 
#' Bourjea, J., Bowen, B.W., Dueñas, R.B., Casale, P., Choudhury, B.C., Costa, A., 
#' Dutton, P.H., Fallabrino, A., Girard, A., Girondot, M., Godfrey, M.H., Hamann, M., 
#' López-Mendilaharsu, M., Marcovaldi, M.A., Mortimer, J.A., Musick, J.A., Nel, R., 
#' Seminoff, J.A., Troëng, S., Witherington, B., Mast, R.B., 2010. Regional 
#' management units for marine turtles: a novel framework for prioritizing 
#' conservation and research across multiple scales. PLoS One 5, e15465.\cr
#' Variance for each value is additive based on both the observed SE (in the RMU.data 
#' object) and a constant value dependent on the rookery when model.SD is equal to 
#' "Rookery-constant". The value is a global constant when model.SD is "global-constant". 
#' The value is proportional to the observed number of nests when model.SD is 
#' "global-proportional" with aSD_*observed+SD_ with aSD_ and SD_ being fitted 
#' values. This value is fixed to zero when model.SD is "Zero".\cr
#' If replicate.CI is 0, no CI is estimated, and only point estimation is returned.
#' @examples
#' \dontrun{
#' library("phenology")
#' RMU.names.AtlanticW <- data.frame(mean=c("Yalimapo.French.Guiana", 
#'                                          "Galibi.Suriname", 
#'                                          "Irakumpapy.French.Guiana"), 
#'                                  se=c("se_Yalimapo.French.Guiana", 
#'                                       "se_Galibi.Suriname", 
#'                                       "se_Irakumpapy.French.Guiana"), 
#'                                  density=c("density_Yalimapo.French.Guiana", 
#'                                            "density_Galibi.Suriname", 
#'                                            "density_Irakumpapy.French.Guiana"), 
#'                                            stringsAsFactors = FALSE)
#' data.AtlanticW <- data.frame(Year=c(1990:2000), 
#'       Yalimapo.French.Guiana=c(2076, 2765, 2890, 2678, NA, 
#'                                6542, 5678, 1243, NA, 1566, 1566),
#'       se_Yalimapo.French.Guiana=c(123.2, 27.7, 62.5, 126, NA, 
#'                                  230, 129, 167, NA, 145, 20),
#'       density_Yalimapo.French.Guiana=rep("dnorm", 11), 
#'       Galibi.Suriname=c(276, 275, 290, NA, 267, 
#'                        542, 678, NA, 243, 156, 123),
#'       se_Galibi.Suriname=c(22.3, 34.2, 23.2, NA, 23.2, 
#'                            4.3, 2.3, NA, 10.3, 10.1, 8.9),
#'       density_Galibi.Suriname=rep("dnorm", 11), 
#'       Irakumpapy.French.Guiana=c(1076, 1765, 1390, 1678, NA, 
#'                                3542, 2678, 243, NA, 566, 566),
#'       se_Irakumpapy.French.Guiana=c(23.2, 29.7, 22.5, 226, NA, 
#'                                  130, 29, 67, NA, 15, 20), 
#'       density_Irakumpapy.French.Guiana=rep("dnorm", 11), stringsAsFactors = FALSE
#'       )
#'                            
#' cst <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'                colname.year="Year", model.trend="Constant", 
#'                model.SD="Zero")
#'                
#' out.CI.Cst <- CI.RMU(result=cst)
#' 
#' 
#'
#' cst <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'                colname.year="Year", model.trend="Constant", 
#'                model.SD="Zero", 
#'                control=list(trace=1, REPORT=100, maxit=500, parscale = c(3000, -0.2, 0.6)))
#'                
#' cst <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'                colname.year="Year", model.trend="Constant", 
#'                model.SD="Zero", method=c("Nelder-Mead","BFGS"), 
#'                control = list(trace = 0, REPORT = 100, maxit = 500, 
#'                parscale = c(3000, -0.2, 0.6)))
#' expo <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'                colname.year="Year", model.trend="Exponential", 
#'                model.SD="Zero", method=c("Nelder-Mead","BFGS"), 
#'                control = list(trace = 0, REPORT = 100, maxit = 500, 
#'                parscale = c(6000, -0.05, -0.25, 0.6)))
#' YS <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'              colname.year="Year", model.trend="Year-specific", method=c("Nelder-Mead","BFGS"), 
#'              model.SD="Zero")
#' YS1 <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'              colname.year="Year", model.trend="Year-specific", method=c("Nelder-Mead","BFGS"), 
#'              model.SD="Zero", model.rookeries="First-order")
#' YS1_cst <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'              colname.year="Year", model.trend="Year-specific", 
#'              model.SD="Constant", model.rookeries="First-order", 
#'              parameters=YS1$par, method=c("Nelder-Mead","BFGS"))
#' YS2 <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'              colname.year="Year", model.trend="Year-specific",
#'              model.SD="Zero", model.rookeries="Second-order", 
#'              parameters=YS1$par, method=c("Nelder-Mead","BFGS"))
#' YS2_cst <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'              colname.year="Year", model.trend="Year-specific",
#'              model.SD="Constant", model.rookeries="Second-order", 
#'              parameters=YS1_cst$par, method=c("Nelder-Mead","BFGS"))
#'                
#' compare_AIC(Constant=cst, 
#'             Exponential=expo, 
#'             YearSpecific=YS)
#' 
#' compare_AIC(YearSpecific_ProportionsFirstOrder_Zero=YS1,
#' YearSpecific_ProportionsFirstOrder_Constant=YS1_cst)
#' 
#' compare_AIC(YearSpecific_ProportionsConstant=YS,
#'            YearSpecific_ProportionsFirstOrder=YS1,
#'            YearSpecific_ProportionsSecondOrder=YS2)
#'            
#' compare_AIC(YearSpecific_ProportionsFirstOrder=YS1_cst,
#'            YearSpecific_ProportionsSecondOrder=YS2_cst)
#' 
#' plot(cst, main="Use of different beaches along the time", what="total")
#' plot(expo, main="Use of different beaches along the time", what="total")
#' plot(YS2_cst, main="Use of different beaches along the time", what="total")
#' 
#' plot(YS1, main="Use of different beaches along the time")
#' plot(YS1_cst, main="Use of different beaches along the time")
#' plot(YS1_cst, main="Use of different beaches along the time", what="numbers")
#' 
#' # Gamma distribution should be used for MCMC outputs
#' 
#' RMU.names.AtlanticW <- data.frame(mean=c("Yalimapo.French.Guiana", 
#'                                          "Galibi.Suriname", 
#'                                          "Irakumpapy.French.Guiana"), 
#'                                  se=c("se_Yalimapo.French.Guiana", 
#'                                       "se_Galibi.Suriname", 
#'                                       "se_Irakumpapy.French.Guiana"), 
#'                                  density=c("density_Yalimapo.French.Guiana", 
#'                                            "density_Galibi.Suriname", 
#'                                            "density_Irakumpapy.French.Guiana"))
#'                                            
#' data.AtlanticW <- data.frame(Year=c(1990:2000), 
#'       Yalimapo.French.Guiana=c(2076, 2765, 2890, 2678, NA, 
#'                                6542, 5678, 1243, NA, 1566, 1566),
#'       se_Yalimapo.French.Guiana=c(123.2, 27.7, 62.5, 126, NA, 
#'                                  230, 129, 167, NA, 145, 20),
#'       density_Yalimapo.French.Guiana=rep("dgamma", 11), 
#'       Galibi.Suriname=c(276, 275, 290, NA, 267, 
#'                        542, 678, NA, 243, 156, 123),
#'       se_Galibi.Suriname=c(22.3, 34.2, 23.2, NA, 23.2, 
#'                            4.3, 2.3, NA, 10.3, 10.1, 8.9),
#'       density_Galibi.Suriname=rep("dgamma", 11), 
#'       Irakumpapy.French.Guiana=c(1076, 1765, 1390, 1678, NA, 
#'                                3542, 2678, 243, NA, 566, 566),
#'       se_Irakumpapy.French.Guiana=c(23.2, 29.7, 22.5, 226, NA, 
#'                                  130, 29, 67, NA, 15, 20), 
#'       density_Irakumpapy.French.Guiana=rep("dgamma", 11)
#'       )
#' cst <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'                colname.year="Year", model.trend="Constant", 
#'                model.SD="Zero")
#' }
#' @export

CI.RMU <- function(result=stop("A result obtained from fitRMU is necessary"), 
                   resultMCMC=NULL, 
                   chain=1, 
                   replicate.CI=10000, 
                   regularThin = TRUE, 
                   silent=FALSE) {
  
  #  result=NULL; resultMCMC=NULL; chain=1; replicate.CI=10000; silent=FALSE; regularThin = TRUE
  
  result$RMU.names$mean <- as.character(result$RMU.names$mean)
  rownames(result$RMU.names) <- result$RMU.names$mean
  if (any(colnames(result$RMU.names) == "se")) {
    result$RMU.names$se <- as.character(result$RMU.names$se)
  }
  if (any(colnames(result$RMU.names) == "density")) {
    result$RMU.names$density <- as.character(result$RMU.names$density)
  }
  
  df_random <- NULL
  hessian <- result$hessian
  pfixed <- c(result$fixed.parameters.initial, result$fixed.parameters.computing)
  totpar <- c(result$par, pfixed)
                                       
  # J'ai une matrice hessienne et pas de MCMC et replicate.CI ne vaut pas 0
  if (((!is.null(hessian)) & is.null(resultMCMC)) & (replicate.CI != 0)) {
    sigma <- try(solve(hessian), silent = TRUE)
    if (!inherits(sigma, "try-error")) {
      s. <- svd(sigma)
      R <- t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
      df_random <- matrix(rnorm(replicate.CI * ncol(sigma)), nrow = replicate.CI, byrow = TRUE) %*% R
      df_random <- sweep(df_random, 2, result$par[rownames(hessian)], "+")
      colnames(df_random) <- rownames(hessian)
    } else {
      # J'ai une erreur sur l'inversion de la matrice hessienne; 
      # Je prends les SE
      
      se <- result$SE
      
      df_random <- matrix(data = NA, ncol=length(se), nrow=replicate.CI)
      colnames(df_random) <- names(se)
      
      for (i in names(se)) {
        df_random[, i] <- rnorm(replicate.CI, result$par[i], se[i])
      }
      
    }
    
    if (!is.null(pfixed)) {
      ajouter <- matrix(rep(pfixed, replicate.CI), 
                        nrow=replicate.CI, byrow=TRUE,
                        dimnames = list(c(NULL), names(pfixed)))
      df_random <- cbind(df_random, ajouter)
    }
  }
  
  if ((!is.null(resultMCMC)) & (replicate.CI != 0)) {
    if (replicate.CI == "all") {
      replicate.CI <- nrow(resultMCMC$resultMCMC[[chain]])
    }
    repl <- FALSE
    if (replicate.CI > nrow(resultMCMC$resultMCMC[[chain]])) repl <- TRUE
    
    df_random <- resultMCMC$resultMCMC[[chain]][sample(1:nrow(resultMCMC$resultMCMC[[chain]]), replicate.CI, replace = repl), ]
    if (!is.null(pfixed)) {
      ajouter <- matrix(rep(pfixed, replicate.CI), 
                        nrow=replicate.CI, byrow=TRUE,
                        dimnames = list(c(NULL), names(pfixed)))
      df_random <- cbind(df_random, ajouter)
    }
  }
  
  if (is.null(df_random)) {
    # Je n'ai pas de moyen de calculer le CI
    if (!silent) warning('No confidence interval is calculated')
    replicate.CI_ec <- 1
    df_random <- matrix(data=totpar, nrow=1, 
                        dimnames = list(c(NULL), names(totpar)))
  } else {
    replicate.CI_ec <- replicate.CI
  }
  
  
  # Mettre aSD_ et SD_ en positif
  
  if (any(substr(colnames(df_random), 1, 4) == "aSD_")) {
    df_random[, substr(colnames(df_random), 1, 4) == "aSD_"] <- ifelse(df_random[, substr(colnames(df_random), 1, 4) == "aSD_"] < 0, 0, df_random[, substr(colnames(df_random), 1, 4) == "aSD_"])
  }
  
  if (any(substr(colnames(df_random), 1, 3) == "SD_")) {
    df_random[, substr(colnames(df_random), 1, 3) == "SD_"] <- ifelse(df_random[, substr(colnames(df_random), 1, 3) == "SD_"] < 0, 0, df_random[, substr(colnames(df_random), 1, 3) == "SD_"])
  }
  
  if (any(substr(colnames(df_random), 1, 2) == "T_")) {
    df_random[, substr(colnames(df_random), 1, 2) == "T_"] <- ifelse(df_random[, substr(colnames(df_random), 1, 2) == "T_"] < 0, 0, df_random[, substr(colnames(df_random), 1, 2) == "T_"])
  }
  
  if (any((substr(colnames(df_random), 1, 1) == "a") & (substr(colnames(df_random), 3, 1) == "_"))) {
    df_random[, (substr(colnames(df_random), 1, 1) == "a") & (substr(colnames(df_random), 3, 1) == "_")] <- ifelse(df_random[, (substr(colnames(df_random), 1, 1) == "a") & (substr(colnames(df_random), 3, 1) == "_")] < 0, 0, df_random[, (substr(colnames(df_random), 1, 1) == "a") & (substr(colnames(df_random), 3, 1) == "_")])
  }
  
  
  
  nbeach <- nrow(result$RMU.names)
  nabeach <- as.character(result$RMU.names[, "mean"])
  
  years <- as.character(min(as.numeric(result$RMU.data[, result$colname.year])):max(result$RMU.data[, result$colname.year]))
  nyear <- length(years)
  
  
  
  # Dans df_random, j'ai les paramètres
  # Maintenant je vais faire un array avec (replicat, année, plage)
  
  map_prop <- array(data = rep(NA, replicate.CI_ec*nbeach*nyear), 
                    dim=c(replicate.CI_ec, nyear, nbeach), 
                    dimnames = list(c(NULL), years, nabeach))
  
  map_number <- map_prop
  
  cumulTot <- NULL
  
  errmissing <- FALSE
  
  for (rep in 1:replicate.CI_ec) {
    
    # D'abord je génère le modèle des proportions par site
    x <- df_random[rep, ]
    
    
    La0 <- x[paste0("a0_", nabeach)]
    La1 <- x[paste0("a1_", nabeach)]
    La2 <- x[paste0("a2_", nabeach)]
    names(La2) <- names(La1) <- names(La0) <- nabeach
    # Dans map, j'ai une matrice avec les plages en colonnes et les années en ligne
    
    # Nombre d'années au total, pas seulement ceux où on a des données
    for (beach in nabeach) {
      map_prop[rep, , beach] <- abs(La2[beach]) * (1:nyear)^2 + abs(La1[beach]) * (1:nyear) + abs(La0[beach])
    }
    
    for (j in 1:nyear) {
      map_prop[rep, j, ] <- map_prop[rep, j, ] / sum(map_prop[rep, j, ])
    }
    
    map_number[rep, , ] <- map_prop[rep, , ]
    
    # J'ai dans map_prop les proportions
    # Maintenant je mets le nombre
    
    
    if (result$model.trend == "year-specific") {
      Tot <- rep(NA, nyear)
      names(Tot) <- years
      Tot[years] <- abs(x[paste0("T_", years)])
      if (any(is.na(Tot))) {
        if (!silent) errmissing <- TRUE
      }
    }
    if (result$model.trend == "constant") {
      Tot <- rep(abs(x["T_"]), nyear)
    }
    if (result$model.trend == "exponential") {
      Tot <- abs(x["T_"]) * exp(x["r"] * (1:nyear))
    }
    
    for (j in 1:nyear) {
      map_number[rep, j, ] <- map_number[rep, j, ] * Tot[j]
    }
    
    cumulTot <- c(cumulTot, Tot)
  }
  
  if (errmissing) warning("Some years are missing; the option year-specific cannot be used safely")
  
  
  tot2 <- matrix(cumulTot, nrow=nyear, byrow = FALSE)
  
  dfTot <- apply(tot2, MARGIN=1, FUN=function(x) quantile(x, probs=c(0.025, 0.5, 0.975), na.rm = TRUE))
  dfTot <- rbind(dfTot, Mean=c(apply(tot2, MARGIN=1, FUN=mean)))
  dfTot <- rbind(dfTot, SD=c(apply(tot2, MARGIN=1, FUN=sd)))
  colnames(dfTot) <- years
  
  map_prop_synthesis <- array(data = rep(NA, 5*nbeach*nyear), 
                              dim=c(5, nyear, nbeach), 
                              dimnames = list(c("2.5%", "50%", "97.5%", "Mean", "SD"), years, nabeach))
  map_number_synthesis <- map_prop_synthesis
  
  for (beach in nabeach) {
    for (y in years) {
      n <- map_prop[, y, beach]
      nq <- quantile(n, probs=c(0.025, 0.5, 0.975), na.rm = TRUE)
      nm <- mean(n, na.rm = TRUE)
      ns <- sd(n, na.rm = TRUE)
      map_prop_synthesis[, y, beach] <- c(nq, nm, ns)
      n <- map_number[, y, beach]
      nq <- quantile(n, probs=c(0.025, 0.5, 0.975), na.rm = TRUE)
      nm <- mean(n, na.rm = TRUE)
      ns <- sd(n, na.rm = TRUE)
      map_number_synthesis[, y, beach] <- c(nq, nm, ns)
    }
  }
  
  # Les observations sont dans result$RMU.data
  # Les noms de colonne sont dans result$RMU.names avec mean, se, et density
  
  map_number_both <- map_number
  
  for (beach in nabeach) {
    for (y in years) {
      
      if ((y %in% as.character(result$RMU.data[, result$colname.year])) & (beach %in% colnames(result$RMU.data))) {
        pos_y <- which(result$RMU.data[, result$colname.year] == y)
        mean_number <- result$RMU.data[pos_y, beach]
      } else {
        mean_number <- NA
      }
      
      if (!is.na(mean_number)) {
        # J'ai une valeur
        
        if (replicate.CI == 0) {
          v <- mean_number
        } else {
          
          if (any(colnames(result$RMU.names) == "density")) {
            density <- result$RMU.data[pos_y, result$RMU.names[beach, "density"]]
          } else {
            density <- ""
          }
          
          if (any(colnames(result$RMU.names) == "se")) {
            se <- result$RMU.data[pos_y, result$RMU.names[beach, "se"]]
          } else {
            se <- 0
          }
          
          if (density == "dnorm") {
            v <- rnorm(replicate.CI, mean=mean_number, sd=se)
          }
          if ((density == "dgamma") & (se != 0)) {
            scale <- se^2/mean_number
            shape <- mean_number*mean_number/(se^2)
            rate <- 1/scale
            v <- rgamma(replicate.CI, shape=shape, rate=rate)
          }
          if ((density == "")  | (se == 0)) {
            v <- rep(mean_number, replicate.CI)
          }
          
          # print(y);print(b); print(mean(map_number_both[, y, b])); print(mean(v))
        }
        # Correction 23/6/2021; c'est y et non pos_y
        map_number_both[, y, beach] <- v
        
      }
    }
  }
  
  # Je calcule maintenant les stats totales
  cumulTot <- NULL
  for(y in years) {
    dfi <- map_number_both[, y,  , drop = FALSE]
    cumulTot <- c(cumulTot, rowSums(x=dfi, dims = 1, na.rm = TRUE))
  }
  
  tot <- matrix(cumulTot, nrow=nyear, byrow = TRUE)
  
  dfTot_both <- apply(tot, MARGIN=1, FUN=function(x) quantile(x, probs=c(0.025, 0.5, 0.975), na.rm = TRUE))
  dfTot_both <- rbind(dfTot_both, Mean=c(apply(tot, MARGIN=1, FUN=mean, na.rm=TRUE)))
  dfTot_both <- rbind(dfTot_both, SD=c(apply(tot, MARGIN=1, FUN=sd, na.rm=TRUE)))
  colnames(dfTot_both) <- years
  
  # Je vérifie que si j'ai des années manquante en year-specific, ça doit être NA
  if (result$model.trend=="year-specific") {
    dfTot_both[, !(colnames(dfTot_both) %in% years)] <- NA
  }
  
  map_number_synthesis_both <- array(data = rep(NA, 5*nbeach*nyear), 
                                     dim=c(5, nyear, nbeach), 
                                     dimnames = list(c("2.5%", "50%", "97.5%", "Mean", "SD"), years, nabeach))
  
  for (beach in nabeach) {
    for (y in years) {
      n <- map_number_both[, y, beach,drop = FALSE]
      nq <- quantile(n, probs=c(0.025, 0.5, 0.975), na.rm = TRUE)
      nm <- mean(n, na.rm = TRUE)
      ns <- sd(n, na.rm = TRUE)
      map_number_synthesis_both[, y, beach] <- c(nq, nm, ns)
    }
  }
  
  return(list(Total=dfTot, 
              Total_both=dfTot_both, 
              Proportions=map_prop_synthesis, 
              Numbers=map_number_synthesis, 
              Numbers_both=map_number_synthesis_both))
}
