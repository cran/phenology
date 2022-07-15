#' fitRMU is used to estimate missing information when several linked values are observed along a timeseries
#' @title Adjust incomplete timeseries with various constraints.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Return a list with the results from optim and synthesis for proportions and numbers
#' @param RMU.data A data.frame with a column Year (the name is defined in colname.year) and one to three columns per rookery defined in RMU.names
#' @param years.byrow If TRUE, the RMU.data data.frame is organized with years in rows
#' @param model.trend Can be Constant, Exponential or Year-specific
#' @param model.rookeries Description temporal change in rookeries proportion. It be Constant, First-order or Second-order
#' @param model.SD Can be Zero, Global-constant, Global-proportional or Rookery-constant. See description.
#' @param parameters  Parameters to fit
#' @param fixed.parameters Parameters that are fixed
#' @param hessian If TRUE, the hessian matrix is calculated and then the standard error of parameters.
#' @param SE Parameters SE for example from fitRMU_MHmcmc()
#' @param replicate.CI Number of replicates to estimate CI of proportion for each rookery
#' @param colname.year Name of the column to be used as time index
#' @param RMU.names A dataframe with one to three columns indicating name of columns for mean, standard deviation, and distribution for roockeris
#' @param method Methods to be used by optim()
#' @param control List of controls for optim()
#' @param itnmax A vector with maximum iterations for each method.
#' @param cptmax.optim How many times optim can be ran when likelihood is better.
#' @param limit.cpt.optim Limit to consider that likelihood is better.
#' @param maxL If an error is produced during the estimation of likelihood, replace -Ln L by this value
#' @family Fill gaps in RMU
#' @description The data must be a data.frame with the first column being years 
#' and two columns for each beach: the average and the se for the estimate.\cr
#' The correspondence between mean, se and density for each rookery are given in the RMU.names data.frame.\cr
#' This data.frame must have a column named mean, another named se and a third named density. If 
#' no sd column exists, no sd will be considered for the series and if no density column exists, it 
#' will be considered as being "dnorm" (Gaussian distribution).\cr
#' The aggregated number of nests and its confidence interval can be obtained using CI.RMU().\cr
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
#' object) and a value.\cr
#' The value is a global constant when model.SD is "global-constant". \cr
#' The value is proportional to the observed number of nests when model.SD is 
#' "global-proportional" with aSD_*observed+SD_ with aSD_ and SD_ being fitted 
#' values. This value is fixed to zero when model.SD is "Zero".\cr
#' The value is dependent on the rookery when model.SD is equal to 
#' "Rookery-constant" or "Rookery-proportional" with a similar formula as previously 
#' described for "global".\cr
#' if method is NULL, it will simply return the names of required parameters.
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
#'                                            "density_Irakumpapy.French.Guiana"))
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
#'       density_Irakumpapy.French.Guiana=rep("dnorm", 11)
#'       )
#'                            
#' cst <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'                colname.year="Year", model.trend="Constant", 
#'                model.SD="Zero")
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
#' compare_AIC(Constant=cst, Exponential=expo, 
#' YearSpecific=YS)
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
#' # Example of different types of plots
#' plot(cst, main="Use of different beaches along the time", what="total", 
#'      ylim=c(0, 4000))
#' plot(cst, main="Use of different beaches along the time", what = "proportions", 
#'      replicate.CI=0)
#' plot(cst, main="Use of different beaches along the time", what = "numbers", 
#'      aggregate="model", ylim=c(0, 4000), replicate.CI=0)
#' plot(cst, main="Use of different beaches along the time", what = "numbers", 
#'      aggregate="both", ylim=c(0, 11000), replicate.CI=0)
#'      
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
#'                                            "density_Irakumpapy.French.Guiana"), 
#'                                            stringsAsFactors = FALSE)
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
#'       density_Irakumpapy.French.Guiana=rep("dgamma", 11), stringsAsFactors = FALSE
#'       )
#' cst <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'                colname.year="Year", model.trend="Constant", 
#'                model.SD="Zero")
#' }
#' @export


#################################################################################
# Ajustement
#################################################################################

fitRMU <- function (RMU.data = stop("data parameter must be provided"), 
                    years.byrow = TRUE, RMU.names = NULL, model.trend = "Constant", 
                    model.rookeries = "Constant", model.SD = "Global-constant", 
                    parameters = NULL, fixed.parameters = NULL, SE = NULL, 
                    method = c("Nelder-Mead", "BFGS"), 
                    control = list(trace = 1), 
                    itnmax = c(1500, 1500), cptmax.optim=100, limit.cpt.optim=1E-5, 
                    hessian = TRUE, replicate.CI = 1000, 
                    colname.year = "Year", maxL = 1e+09) 
{
  
  #  RMU.data = NULL; years.byrow = TRUE; RMU.names = NULL
  #  model.trend = "Constant"
  #  model.rookeries = "Constant"; model.SD = "Global-constant"; cptmax.optim=100
  #  parameters = NULL; fixed.parameters = NULL; SE = NULL
  #  method = c("Nelder-Mead", "BFGS")
  #  control = list(trace = 1, REPORT = 100, maxit = 1500)
  #  itnmax = c(1500, 1500); hessian = TRUE
  #  replicate.CI = 1000; colname.year = "Year"; maxL = 1e+09
  
  if (!years.byrow) RMU.data <- t(RMU.data)
  
  fpar <- fixed.parameters
  model.trend <- tolower(model.trend)
  model.rookeries <- tolower(model.rookeries)
  model.SD <- tolower(model.SD)
  model.trend <- match.arg(model.trend, choices = c("year-specific", 
                                                    "constant", "exponential"))
  model.rookeries <- match.arg(model.rookeries, choices = c("first-order", 
                                                            "constant", "second-order"))
  model.SD <- match.arg(model.SD, choices = c("global-constant", "rookery-proportional", 
                                              "global-proportional", "rookery-constant", "zero"))
  nm <- colnames(RMU.data)
  index.year <- which(nm == colname.year)
  
  if (all(colnames(RMU.names) != "density")) {
    index.density <- NULL
  } else {
    index.density <- match(RMU.names$density, nm)
  }
  
  index.mean <- match(RMU.names$mean, nm)
  if (all(colnames(RMU.names) != "se")) {
    index.se <- NULL
  }    else {
    index.se <- match(RMU.names$se, nm)
  }
  d <- RMU.data[, index.mean]
  nabeach <- colnames(d)
  nbeach <- length(nabeach)
  nyear <- dim(RMU.data)[1]
  nayear <- RMU.data[, index.year]
  
  print(paste0("Number of available data: ", sum(!is.na(d))))
  index <- list(year = index.year, mean = index.mean, se = index.se, 
                density = index.density, 
                colnames = nabeach, nyear = nyear, nbeach = nbeach, maxL = maxL)
  if (any(substr(nabeach, 1, 2) == "T_") | any(substr(nabeach, 
                                                      1, 3) == "a0_") | any(substr(nabeach, 1, 3) == "a1_") | 
      any(substr(nabeach, 1, 3) == "a2_") | any(substr(nabeach, 
                                                       1, 3) == "SD_") | any(nabeach == "r")) {
    stop("Names of rookeries cannot begin with T_, a0_, a1_, a2_ or SD_ and cannot be r. Please change them.")
  }
  
  if (!is.null(index.se) & !is.null(index.density)) {
    if (!all(c(index.year, index.mean, index.se, index.density) %in% seq_along(nm))) 
      stop("check the correspondance between names of columns and RMU.names or colname.year")
  }
  
  if (!is.null(index.se)) {
    if (!all(c(index.year, index.mean, index.se) %in% seq_along(nm))) 
      stop("check the correspondance between names of columns and RMU.names or colname.year")
  }
  
  if (!is.null(index.density)) {
    if (!all(c(index.year, index.mean, index.density) %in% seq_along(nm))) 
      stop("check the correspondance between names of columns and RMU.names or colname.year")
  }
  
  if (!all(c(index.year, index.mean) %in% seq_along(nm))) 
    stop("check the correspondance between names of columns and RMU.names or colname.year")
  
  
  Tot <- NULL
  LikelihoodRMU <- getFromNamespace(".LikelihoodRMU", ns = "phenology")
  # inv.p.multinomial <- getFromNamespace(".inv.p.multinomial", 
  #                                      ns = "phenology")
  # p.multinomial <- getFromNamespace(".p.multinomial", ns = "phenology")
  if (model.trend == "year-specific") {
    intT <- d
    cm <- colMeans(d, na.rm = TRUE)
    for (beach in 1:nbeach) intT[is.na(intT[, beach]), beach] <- cm[beach]
    Tot <- rowSums(intT)
    names(Tot) <- paste0("T_", RMU.data[, index.year])
    SD <- apply(intT, 2, function(x) sd(x, na.rm = TRUE))/8
    names(SD) <- paste0("SD_", names(SD))
  }
  if (model.trend == "constant") {
    Tot <- c(T_ = sum(colMeans(d, na.rm = TRUE), na.rm = TRUE))
    SD <- apply(d, 2, function(x) sd(x, na.rm = TRUE))/2
    names(SD) <- paste0("SD_", names(SD))
  }
  if (model.trend == "exponential") {
    Tot <- c(T_ = sum(colMeans(d, na.rm = TRUE), na.rm = TRUE), 
             r = 0.001)
    SD <- apply(d, 2, function(x) sd(x, na.rm = TRUE))/2
    names(SD) <- paste0("SD_", names(SD))
  }
  SD[is.na(SD)] <- 1
  SD[SD == 0] <- 1
  
  #26/10/2019
  aSD <- SD
  names(aSD) <- paste0("a", names(SD))
  
  p <- colMeans(d, na.rm = TRUE)
  p[p == 0] <- 1E-5
  p <- p/p[1]
  
  a0 <- p
  a1 <- a0
  a2 <- a0
  a2[] <- 0
  a1[] <- 0
  names(a0) <- paste0("a0_", names(p))
  names(a1) <- paste0("a1_", names(p))
  names(a2) <- paste0("a2_", names(p))
  if (model.rookeries == "constant") {
    x <- c(Tot, a0[-1])
    fixed.parameters <- c(fixed.parameters, a1, a2, a0[1])
  }    else {
    if (model.rookeries == "first-order") {
      x <- c(Tot, a0[-1], a1[-1])
      fixed.parameters <- c(fixed.parameters, a2, a0[1], a1[1])
    }        else {
      x <- c(Tot, a0[-1], a1[-1], a2[-1])
      fixed.parameters <- c(fixed.parameters, a2[1], a0[1], a1[1])
    }
  }
  if (model.SD == "zero") {
    SD[] <- 0
    fixed.parameters <- c(fixed.parameters, SD)
  }
  if ((model.SD == "global-constant") | (model.SD == "global-proportional")) {
    x <- c(x, SD_ = max(d, na.rm = TRUE)/10)
    if ((model.SD == "global-proportional")) {
      x <- c(x, aSD_ = 0.01)
    }
  }
  if ((model.SD == "rookery-constant")  | (model.SD == "rookery-proportional")) {
    x <- c(x, SD)
    if ((model.SD == "rookery-proportional")) {
      x <- c(x, aSD)
    }
  }
  if (any(!is.null(parameters))) {
    for (i in 1:length(parameters)) {
      x <- x[names(x) != names(parameters[i])]
    }
    x <- c(x, parameters)
  }
  if (!is.null(fixed.parameters)) {
    for (i in 1:length(fixed.parameters)) {
      x <- x[names(x) != names(fixed.parameters[i])]
    }
  }
  if (length(x) >= sum(!is.na(d))) {
    stop("More parameters to be fitted than data. Use MCMC rather.")
  }
  
  Lprecedent <- +Inf
  cpt.optim <- 1
  
  # print(d(x))
  
  if (is.null(method)) return(list(fitted.parameters=x, fixed.parameters=fixed.parameters))
  
  nm <- names(x)
  
  repeat {
    # rep(1, length(x))
    scale.factor <- x
    scale.factor <- ifelse(scale.factor == 0, 1, scale.factor)
    
    itnmax <- rep(itnmax, length(method))[1:length(method)]
    
    for (i in 1:length(method)) {
      result <- optim(par = x, 
                      fn = LikelihoodRMU, gr = NULL, fixed.parameters = fixed.parameters, 
                      RMU.data = RMU.data, index = index, model.trend = model.trend, 
                      colname.year=colname.year, RMU.names = RMU.names , 
                      model.SD = model.SD, 
                      method = method[i], 
                      control = modifyList(control, list(maxit = itnmax[i], 
                                                         parscale=scale.factor)), 
                      hessian = FALSE)
      
      x <- result$par
      # C'est utile à cause de la méthode Brent
      names(x) <- nm
    }
    
    # if (any(class(result) == "try-error")) 
    #   stop("An error occurred during the fit. Check initial conditions.")
    # minL <- dim(result)[1]
    # nm <- names(x)
    # x <- result[minL, 1:length(nm), drop = TRUE]
    # x <- unlist(x)
    # names(x) <- nm
    # conv <- result[minL, "convcode"]
    value <- result$value
    x[substr(names(x), 1, 2) == "T_"] <- abs(x[substr(names(x), 
                                                      1, 2) == "T_"])
    x[substr(names(x), 1, 3) == "SD_"] <- abs(x[substr(names(x), 
                                                       1, 3) == "SD_"])
    x[substr(names(x), 1, 4) == "aSD_"] <- abs(x[substr(names(x), 
                                                        1, 4) == "aSD_"])
    x[substr(names(x), 1, 3) == "a0_"] <- abs(x[substr(names(x), 
                                                       1, 3) == "a0_"])
    x[substr(names(x), 1, 3) == "a1_"] <- abs(x[substr(names(x), 
                                                       1, 3) == "a1_"])
    x[substr(names(x), 1, 3) == "a2_"] <- abs(x[substr(names(x), 
                                                       1, 3) == "a2_"])
    if ((abs(Lprecedent-value) < limit.cpt.optim) | (cpt.optim >= cptmax.optim))  break
    Lprecedent <- value
    cpt.optim <- cpt.optim + 1
    print(paste("Convergence is not achieved -LnL=", 
                format(value, digit = floor(log(abs(value))/log(10)) + 3)))
  }
  if (hessian) {
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      stop("numDeriv package is absent; Please install it first")
    }
    message("Estimation of the standard error of parameters. Be patient please.")
    mathessian <- try(getFromNamespace("hessian", ns = "numDeriv")(func = LikelihoodRMU, 
                                                                   fixed.parameters = fixed.parameters, 
                                                                   RMU.data = RMU.data, index = index, 
                                                                   model.trend = model.trend, model.SD = model.SD, x = x, 
                                                                   colname.year=colname.year, RMU.names = RMU.names , 
                                                                   method = "Richardson"), silent = TRUE)
    if (substr(mathessian[1], 1, 5) == "Error") {
      res_se <- rep(NA, length(x))
      names(res_se) <- names(x)
    } else {
      rownames(mathessian) <- colnames(mathessian) <- names(x)
      res_se <- SEfromHessian(mathessian)
    }
  } else {
    warning("Standard errors are not estimated.")
    mathessian <- NULL
    res_se <- rep(NA, length(x))
    names(res_se) <- names(x)
  }
  
  result$hessian <- mathessian
  result$SE <- modifyVector(res_se, SE)
  if (value > 0) {
    print(paste("Convergence is achieved. -LnL=", 
                format(value, digit = floor(log(abs(value))/log(10)) + 3)))
  }    else {
    print(paste("Convergence is achieved. -LnL=", value))
  }
  
  result$AIC <- 2 * value + 2 * length(x)
  print(paste("Parameters=", length(x)))
  
  if (result$AIC > 0) {
    print(paste("AIC=", format(result$AIC, digit = floor(log(abs(result$AIC))/log(10)) + 3)))
  }    else {
    print(paste("AIC=", result$AIC))
  }
  
  result$model.trend <- model.trend
  result$model.rookeries <- model.rookeries
  result$RMU.data <- RMU.data
  result$model.SD <- model.SD
  result$RMU.names <- RMU.names
  result$fixed.parameters.initial <- fpar
  result$fixed.parameters.computing <- fixed.parameters
  result$replicate.CI <- replicate.CI
  result$colname.year <- colname.year
  result <- addS3Class(result, "fitRMU")
  return(result)
}
