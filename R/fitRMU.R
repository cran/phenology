#' fitRMU is used to estimate missing information when several linked values are observed along a timeseries
#' @title Adjust incomplete timeseries with various constraints.
#' @author Marc Girondot
#' @return Return a list with the results from optim and synthesis for proportions and numbers
#' @param RMU.data A data.frame with a column Year (the name is defined in colname.year) and two columns per rookery defined in RMU.names
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
#' @param RMU.names A dataframe with two columns indicating name of columns for mean and standard error for roockerys
#' @param method Methods to be used by optimx()
#' @param control List of controls for optimx()
#' @param itnmax A vector with maximum iterations for each method.
#' @param maxL If an error is produced during the estimation of likelihood, replace -Ln L by this value
#' @family Fill gaps in RMU
#' @description The data must be a data.frame with the first column being years \cr
#' and two columns for each beach: the average and the se for the estimate.\cr
#' The correspondance between mean and se for each rookery are given in the RMU.names data.frame.\cr
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
#' "beach-constant". The value is a global constant when model.SD is "global-constant". 
#' The value is proportional to the observed number of nests when model.SD is 
#' "global-proportional" with aSD_*observed+SD_ with aSD_ and SD_ being fitted 
#' values. This value is fixed to zero when model.SD is "Zero".
#' @family Fill the gaps for RMU
#' @examples
#' \dontrun{
#' library("phenology")
#' RMU.names.AtlanticW <- data.frame(mean=c("Yalimapo.French.Guiana", 
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
#' cst <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'                colname.year="Year", model.trend="Constant", 
#'                model.SD="Zero")
#' cst <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'                colname.year="Year", model.trend="Constant", 
#'                model.SD="Zero", 
#'                control=list(trace=1, REPORT=100, maxit=500, parscale = c(3000, -0.2, 0.6)))
#'                
#' # Example with optimx
#' require("optimx")
#' cst <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'                colname.year="Year", model.trend="Constant", 
#'                model.SD="Zero", optim="optimx", method=c("Nelder-Mead","BFGS"), 
#'                control = list(trace = 0, REPORT = 100, maxit = 500, 
#'                parscale = c(3000, -0.2, 0.6)))
#' expo <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'                colname.year="Year", model.trend="Exponential", 
#'                model.SD="Zero", optim="optimx", method=c("Nelder-Mead","BFGS"), 
#'                control = list(trace = 0, REPORT = 100, maxit = 500, 
#'                parscale = c(6000, -0.05, -0.25, 0.6)))
#' YS <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'              colname.year="Year", model.trend="Year-specific", method=c("Nelder-Mead","BFGS"), 
#'              optim="optimx", model.SD="Zero")
#' YS1 <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'              colname.year="Year", model.trend="Year-specific", method=c("Nelder-Mead","BFGS"), 
#'              optim="optimx", model.SD="Zero", model.rookeries="First-order")
#' YS1_cst <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'              colname.year="Year", model.trend="Year-specific", 
#'              model.SD="Constant", model.rookeries="First-order", 
#'              optim="optimx", parameters=YS1$par, method=c("Nelder-Mead","BFGS"))
#' YS2 <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'              colname.year="Year", model.trend="Year-specific",
#'              model.SD="Zero", model.rookeries="Second-order", 
#'              optim="optimx", parameters=YS1$par, method=c("Nelder-Mead","BFGS"))
#' YS2_cst <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'              colname.year="Year", model.trend="Year-specific",
#'              model.SD="Constant", model.rookeries="Second-order", 
#'              optim="optimx", parameters=YS1_cst$par, method=c("Nelder-Mead","BFGS"))
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
#' barplot_errbar(YS1_cst$proportions[1, ], y.plus = YS1_cst$proportions.CI.0.95[1, ], 
#' y.minus = YS1_cst$proportions.CI.0.05[1, ], las=1, ylim=c(0, 0.7), 
#' main="Proportion of the different rookeries in the region")
#' 
#' plot(cst, main="Use of different beaches along the time", what="total")
#' plot(expo, main="Use of different beaches along the time", what="total")
#' plot(YS2_cst, main="Use of different beaches along the time", what="total")
#' 
#' plot(YS1, main="Use of different beaches along the time")
#' plot(YS1_cst, main="Use of different beaches along the time")
#' plot(YS1_cst, main="Use of different beaches along the time", what="numbers")
#' 
#' }
#' @export


#################################################################################
# Ajustement
#################################################################################

fitRMU<-function (RMU.data = stop("data parameter must be provided"), 
                  years.byrow = TRUE, RMU.names = NULL, model.trend = "Constant", 
                  model.rookeries = "Constant", model.SD = "Global-constant", 
                  parameters = NULL, fixed.parameters = NULL, SE = NULL, method = c("Nelder-Mead", 
                                                                                    "BFGS"), control = list(trace = 1, REPORT = 100, maxit = 1500), 
                  itnmax = c(1500, 1500), hessian = TRUE, replicate.CI = 1000, 
                  colname.year = "Year", maxL = 1e+09) 
{
  if (!years.byrow) 
    RMU.data <- t(RMU.data)
  model.trend <- tolower(model.trend)
  model.rookeries <- tolower(model.rookeries)
  model.SD <- tolower(model.SD)
  model.trend <- match.arg(model.trend, choices = c("year-specific", 
                                                    "constant", "exponential"))
  model.rookeries <- match.arg(model.rookeries, choices = c("first-order", 
                                                            "constant", "second-order"))
  model.SD <- match.arg(model.SD, choices = c("global-constant", 
                                              "global-proportional", "rookery-constant", "zero"))
  nm <- colnames(RMU.data)
  index.year <- which(nm == colname.year)
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
  print(paste0("Number of available data:", sum(!is.na(d))))
  index <- list(year = index.year, mean = index.mean, se = index.se, 
                colnames = nabeach, nyear = nyear, nbeach = nbeach, maxL = maxL)
  if (any(substr(nabeach, 1, 2) == "T_") | any(substr(nabeach, 
                                                      1, 3) == "a0_") | any(substr(nabeach, 1, 3) == "a1_") | 
      any(substr(nabeach, 1, 3) == "a2_") | any(substr(nabeach, 
                                                       1, 3) == "SD_") | any(nabeach == "r")) {
    stop("Names of rookeries cannot begin with T_, a0_, a1_, a2_ or SD_ and cannot be r. Please change them.")
  }
  if (!is.null(index.se)) {
    if (!all(c(index.year, index.mean, index.se) %in% seq_along(nm))) 
      stop("check the correspondance between names of columns and RMU.names or colname.year")
  }    else {
    if (!all(c(index.year, index.mean) %in% seq_along(nm))) 
      stop("check the correspondance between names of columns and RMU.names or colname.year")
  }
  Tot <- NULL
  LikelihoodRMU <- getFromNamespace(".LikelihoodRMU", ns = "phenology")
  inv.p.multinomial <- getFromNamespace(".inv.p.multinomial", 
                                        ns = "phenology")
  p.multinomial <- getFromNamespace(".p.multinomial", ns = "phenology")
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
             r = 0)
    SD <- apply(d, 2, function(x) sd(x, na.rm = TRUE))/2
    names(SD) <- paste0("SD_", names(SD))
  }
  SD[is.na(SD)] <- 1
  SD[SD == 0] <- 1
  p <- colMeans(d, na.rm = TRUE)
  p <- p/sum(p)
  p <- inv.p.multinomial(p)
  p <- log((1/p) - 1)
  a0 <- p
  a1 <- p
  a2 <- p
  a2[] <- 0
  a1[] <- 0
  names(a0) <- paste0("a0_", names(p))
  names(a1) <- paste0("a1_", names(p))
  names(a2) <- paste0("a2_", names(p))
  if (model.rookeries == "constant") {
    x <- c(Tot, a0)
    fixed.parameters <- c(fixed.parameters, a1, a2)
  }    else {
    if (model.rookeries == "first-order") {
      x <- c(Tot, a0, a1)
      fixed.parameters <- c(fixed.parameters, a2)
    }        else {
      x <- c(Tot, a0, a1, a2)
    }
  }
  if (model.SD == "zero") {
    SD[] <- 0
    fixed.parameters <- c(fixed.parameters, SD)
  }
  if ((model.SD == "global-constant") | (model.SD == "global-proportional")) {
    x <- c(x, SD_ = 1)
    if ((model.SD == "global-proportional")) {
      x <- c(x, aSD_ = 1)
    }
  }
  if (model.SD == "rookery-constant") {
    x <- c(x, SD)
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
  scale.factor <- rep(1, length(x))
  repeat {
    result <- try(suppressWarnings(optimx::optimx(par = x, 
                                                  fn = LikelihoodRMU, gr = NULL, fixed = fixed.parameters, 
                                                  RMU.data = RMU.data, index = index, model.trend = model.trend, 
                                                  model.SD = model.SD, method = method, itnmax = itnmax, 
                                                  control = modifyList(control, list(dowarn = FALSE, 
                                                                                     follow.on = TRUE)), hessian = FALSE)), silent = TRUE)
    if (any(class(result) == "try-error")) 
      stop("An error occurred during the fit. Check initial conditions.")
    minL <- dim(result)[1]
    nm <- names(x)
    x <- result[minL, nm]
    x <- as.numeric(x)
    names(x) <- nm
    conv <- result[minL, "convcode"]
    value <- result[minL, "value"]
    x[substr(names(x), 1, 2) == "T_"] <- abs(x[substr(names(x), 
                                                      1, 2) == "T_"])
    x[substr(names(x), 1, 3) == "SD_"] <- abs(x[substr(names(x), 
                                                       1, 3) == "SD_"])
    if (conv == 0) 
      break
    print(paste("Convergence is not achieved -LnL=", format(value, 
                                                            digit = floor(log(value)/log(10)) + 3)))
  }
  if (hessian) {
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      stop("numDeriv package is absent; Please install it first")
    }
    message("Estimation of the standard error of parameters. Be patient please.")
    mathessian <- try(getFromNamespace("hessian", ns = "numDeriv")(func = LikelihoodRMU, 
                                                                   fixed = fixed.parameters, RMU.data = RMU.data, index = index, 
                                                                   model.trend = model.trend, model.SD = model.SD, x = x, 
                                                                   method = "Richardson"), silent = TRUE)
    if (substr(mathessian[1], 1, 5) == "Error") {
      res_se <- rep(NA, length(x))
      names(res_se) <- names(x)
    }        else {
      rownames(mathessian) <- colnames(mathessian) <- names(x)
      res_se <- SEfromHessian(mathessian)
    }
  }    else {
    warning("Standard errors are not estimated.")
    mathessian <- NULL
    res_se <- rep(NA, length(x))
    names(res_se) <- names(x)
  }
  result_list <- list()
  result_list$par <- x
  result_list$value <- value
  result_list$convergence <- conv
  result_list$hessian <- mathessian
  result_list$SE <- modifyVector(res_se, SE)
  if (value > 0) {
    print(paste("Convergence is achieved. -LnL=", format(value, 
                                                         digit = floor(log(value)/log(10)) + 3)))
  }    else {
    print(paste("Convergence is achieved. -LnL=", value))
  }
  result_list$AIC <- 2 * value + 2 * length(x)
  print(paste("Parameters=", length(x)))
  if (result_list$AIC > 0) {
    print(paste("AIC=", format(result_list$AIC, digit = floor(log(result_list$AIC)/log(10)) + 
                                 3)))
  }    else {
    print(paste("AIC=", result_list$AIC))
  }
  x <- c(x, fixed.parameters)
  sepf <- fixed.parameters
  sepf[] <- 0
  res <- modifyVector(result_list$SE, sepf)
  La0 <- x[paste0("a0_", nabeach[paste0("a0_", nabeach) %in% 
                                   names(x)])]
  La1 <- x[paste0("a1_", nabeach[paste0("a1_", nabeach) %in% 
                                   names(x)])]
  La2 <- x[paste0("a2_", nabeach[paste0("a2_", nabeach) %in% 
                                   names(x)])]
  map <- matrix(rep(NA, nyear * (nbeach - 1)), ncol = nbeach - 
                  1)
  for (beach in 1:(nbeach - 1)) {
    map[, beach] <- 1/(1 + exp(La2[beach] * (1:nyear)^2 + 
                                 La1[beach] * (1:nyear) + La0[beach]))
  }
  mapp <- matrix(rep(NA, nyear * nbeach), ncol = nbeach, dimnames = list(year = nayear, 
                                                                         rookery = nabeach))
  for (j in 1:nyear) {
    mapp[j, ] <- p.multinomial(map[j, ])
  }
  result_list$proportions <- mapp
  if (any(is.na(res))) {
    result_list$proportions.CI.0.025 <- NA
    result_list$proportions.CI.0.975 <- NA
  }    else {
    sea0 <- res[paste0("a0_", nabeach[paste0("a0_", nabeach) %in% 
                                        names(x)])]
    sea1 <- res[paste0("a1_", nabeach[paste0("a1_", nabeach) %in% 
                                        names(x)])]
    sea2 <- res[paste0("a2_", nabeach[paste0("a2_", nabeach) %in% 
                                        names(x)])]
    map.La0 <- matrix(rep(NA, replicate.CI * (nbeach - 1)), 
                      ncol = nbeach - 1, dimnames = list(NULL, names(La0)))
    map.La1 <- matrix(rep(NA, replicate.CI * (nbeach - 1)), 
                      ncol = nbeach - 1, dimnames = list(NULL, names(La1)))
    map.La2 <- matrix(rep(NA, replicate.CI * (nbeach - 1)), 
                      ncol = nbeach - 1, dimnames = list(NULL, names(La2)))
    for (beach in seq_along(La0)) {
      map.La0[, beach] <- rnorm(replicate.CI, La0[beach], 
                                sea0[beach])
      map.La1[, beach] <- rnorm(replicate.CI, La1[beach], 
                                sea1[beach])
      map.La2[, beach] <- rnorm(replicate.CI, La2[beach], 
                                sea2[beach])
    }
    mapparray <- array(data = NA, dim = c(nyear, nbeach, 
                                          replicate.CI), dimnames = list(year = NULL, rookery = nabeach, 
                                                                         replicate = NULL))
    for (replicate in 1:replicate.CI) {
      for (year in 1:nyear) {
        p <- 1/(1 + exp(map.La2[replicate, ] * (year)^2 + 
                          map.La1[replicate, ] * (year) + map.La0[replicate, 
                                                                  ]))
        mapparray[year, , replicate] <- p.multinomial(p)
      }
    }
    mapp.CI.0.025 <- matrix(rep(NA, nyear * nbeach), ncol = nbeach, 
                            dimnames = list(year = nayear, rookery = nabeach))
    mapp.CI.0.975 <- matrix(rep(NA, nyear * nbeach), ncol = nbeach, 
                            dimnames = list(year = nayear, rookery = nabeach))
    for (year in 1:nyear) {
      for (beach in 1:nbeach) {
        mapp.CI.0.025[year, beach] <- quantile(mapparray[year, 
                                                         beach, ], probs = 0.025)
        mapp.CI.0.975[year, beach] <- quantile(mapparray[year, 
                                                         beach, ], probs = 0.975)
      }
    }
    result_list$proportions.CI.0.025 <- mapp.CI.0.025
    result_list$proportions.CI.0.975 <- mapp.CI.0.975
  }
  if (model.trend == "year-specific") {
    Tot <- x[paste0("T_", nayear)]
    dtaL_theo <- matrix(rep(Tot, nbeach), ncol = nbeach, 
                        byrow = FALSE, dimnames = list(year = nayear, rookery = nabeach))
  }
  if (model.trend == "constant") {
    Tot <- x["T_"]
    dtaL_theo <- matrix(rep(Tot, nbeach * nyear), ncol = nbeach, 
                        byrow = FALSE, dimnames = list(year = nayear, rookery = nabeach))
  }
  if (model.trend == "exponential") {
    Tot <- x["T_"] * exp(x["r"] * (1:nyear))
    dtaL_theo <- matrix(rep(Tot, nbeach), ncol = nbeach, 
                        byrow = FALSE, dimnames = list(year = nayear, rookery = nabeach))
  }
  result_list$numbers <- dtaL_theo * mapp
  if (any(is.na(res))) {
    result_list$numbers.CI.0.025 <- NA
    result_list$numbers.CI.0.975 <- NA
  }    else {
    if (model.trend == "year-specific") {
      sTot <- res[paste0("T_", nayear)]
      arraydtaL_theo <- array(data = NA, dim = c(nyear, 
                                                 nbeach, replicate.CI), dimnames = list(year = nayear, 
                                                                                        rookery = nabeach, replicate = NULL))
      for (replicate in 1:replicate.CI) {
        Totec <- sapply(seq_along(Tot), function(x) rnorm(1, 
                                                          Tot[x], sTot[x]))
        arraydtaL_theo[, , replicate] <- matrix(rep(Totec, 
                                                    nbeach), ncol = nbeach, byrow = FALSE)
      }
    }
    if (model.trend == "constant") {
      sTot <- res["T_"]
      arraydtaL_theo <- array(data = NA, dim = c(nyear, 
                                                 nbeach, replicate.CI), dimnames = list(year = nayear, 
                                                                                        rookery = nabeach, replicate = NULL))
      Tot <- rnorm(replicate.CI, Tot, sTot)
      for (replicate in 1:replicate.CI) arraydtaL_theo[, 
                                                       , replicate] <- matrix(rep(Tot[replicate], nbeach * 
                                                                                    nyear), ncol = nbeach, byrow = FALSE)
    }
    if (model.trend == "exponential") {
      arraydtaL_theo <- array(data = NA, dim = c(nyear, 
                                                 nbeach, replicate.CI), dimnames = list(year = nayear, 
                                                                                        rookery = nabeach, replicate = NULL))
      for (replicate in 1:replicate.CI) {
        Totec <- rnorm(1, x["T_"], res["T_"]) * exp(rnorm(1, 
                                                          x["r"], res["r"]) * (1:nyear))
        arraydtaL_theo[, , replicate] <- matrix(rep(Totec, 
                                                    nbeach), ncol = nbeach, byrow = FALSE)
      }
    }
    arraydtaL_theo2 <- arraydtaL_theo * mapparray
    numbers.CI.0.025 <- matrix(rep(NA, nyear * nbeach), ncol = nbeach, 
                               dimnames = list(year = nayear, rookery = nabeach))
    numbers.CI.0.975 <- matrix(rep(NA, nyear * nbeach), ncol = nbeach, 
                               dimnames = list(year = nayear, rookery = nabeach))
    for (year in 1:nyear) {
      for (beach in 1:nbeach) {
        numbers.CI.0.025[year, beach] <- quantile(arraydtaL_theo2[year, 
                                                                  beach, ], probs = 0.025)
        numbers.CI.0.975[year, beach] <- quantile(arraydtaL_theo2[year, 
                                                                  beach, ], probs = 0.975)
      }
    }
    result_list$numbers.CI.0.025 <- numbers.CI.0.025
    result_list$numbers.CI.0.975 <- numbers.CI.0.975
  }
  result_list$model.trend <- model.trend
  result_list$RMU.data <- RMU.data
  result_list$model.SD <- model.SD
  result_list$RMU.names <- RMU.names
  result_list$fixed.parameters <- fixed.parameters
  result_list$replicate.CI <- replicate.CI
  result_list$colname.year <- colname.year
  class(result_list) = "fitRMU"
  return(result_list)
}
