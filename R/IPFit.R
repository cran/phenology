#' IPFit fit a model of Internesting Period for marine turtles.
#' @title Fit a model of Internesting Period for marine turtles.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Return a list of class IP with the fit informations and the fitted model.\cr
#' @param x Initial parameters to be fitted
#' @param fixed.parameters Parameters that are fixed.
#' @param data Data as a vector
#' @param method Method to be used by optimx()
#' @param itnmax A vector with maximum iterations for each method.
#' @param control List of controls for optimx()
#' @param hessian Logical to estimate SE of parameters
#' @param parallel If TRUE, will use parallel computing
#' @param verbose If TRUE, show the parameters for each tested model
#' @param model Can be ML for Maximum likelihood or MH for Metropolis Hastings
#' @param parametersMH The priors. See MHalgoGen
#' @param n.iter See MHalgoGen
#' @param n.chains See MHalgoGen
#' @param n.adapt See MHalgoGen
#' @param thin See MHalgoGen
#' @param trace See MHalgoGen
#' @param adaptive See MHalgoGen
#' @param adaptive.lag See MHalgoGen
#' @param adaptive.fun See MHalgoGen
#' @param intermediate See MHalgoGen
#' @param filename See MHalgoGen
#' @description This function fits a model of internesting period using maximum 
#' likelihood or using Metropolis-Hastings algorithm with Bayesian model.\cr
#' The fit using maximum likelihood is not the best strategy because the objective 
#' function is based on a stochastic model (and then a single set of parameters 
#' does not produce exactly the same output each time). 
#' The use of Metropolis-Hastings algorithm (a Markov chain Monte Carlo method)  
#' should be prefered.
#' @family Model of Internesting Period
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' data <- structure(c(`0` = 0, `1` = 47, `2` = 15, `3` = 6, `4` = 5, `5` = 4, 
#'                     `6` = 2, `7` = 5, `8` = 57, `9` = 203, `10` = 205, `11` = 103, 
#'                     `12` = 35, `13` = 24, `14` = 12, `15` = 10, `16` = 13, `17` = 49, 
#'                     `18` = 86, `19` = 107, `20` = 111, `21` = 73, `22` = 47, `23` = 30, 
#'                     `24` = 19, `25` = 17, `26` = 33, `27` = 48, `28` = 77, `29` = 83, 
#'                     `30` = 65, `31` = 37, `32` = 27, `33` = 23, `34` = 24, `35` = 22, 
#'                     `36` = 41, `37` = 42, `38` = 44, `39` = 33, `40` = 39, `41` = 24, 
#'                     `42` = 18, `43` = 18, `44` = 22, `45` = 22, `46` = 19, `47` = 24, 
#'                     `48` = 28, `49` = 17, `50` = 18, `51` = 19, `52` = 17, `53` = 4, 
#'                     `54` = 12, `55` = 9, `56` = 6, `57` = 11, `58` = 7, `59` = 11, 
#'                     `60` = 12, `61` = 5, `62` = 4, `63` = 6, `64` = 11, `65` = 5, 
#'                     `66` = 6, `67` = 7, `68` = 3, `69` = 2, `70` = 1, `71` = 3, `72` = 2, 
#'                     `73` = 1, `74` = 2, `75` = 0, `76` = 0, `77` = 3, `78` = 1, `79` = 0, 
#'                     `80` = 2, `81` = 0, `82` = 0, `83` = 1), Year = "1994", 
#'                     Species = "Dermochelys coriacea", 
#'                     location = "Yalimapo beach, French Guiana", 
#'                     totalnumber = 2526L, class = "IP")
#'   par(mar=c(4, 4, 1, 1)+0.4)
#'   plot(data, xlim=c(0,100))
#'   text(100, 190, labels=bquote(italic(.(attributes(data)$Species))), pos=2)
#'   text(100, 150, labels=attributes(data)$location, pos=2, cex=0.8)
#'   text(100, 110, labels=paste0(as.character(attributes(data)$totalnumber), " females"), pos=2)
#'   
#' ######### Fit using Maximum-Likelihood
#' 
#' par <- c(meanIP = 9.8229005713237623, 
#'          sdIP = 0.079176011861863474, 
#'          minIP = 6.8128364577569309, 
#'          pAbort = 1.5441529841959203, 
#'          meanAbort = 2.7958742380756121, 
#'          sdAbort = 0.99370406770777175, 
#'          pCapture = -0.80294884905867658, 
#'          meanECF = 4.5253772889275758, 
#'          sdECF = 0.20334743335612529)
#' 
#' fML <- IPFit(x=par, 
#'              fixed.parameters=c(N=20000),
#'              data=data, 
#'              verbose=FALSE, 
#'              model="ML")
#' 
#' # Plot the fitted ECF
#' plot(fML, result="ECF")
#' 
#' # Plot the Internesting Period distribution
#' plot(fML, result="IP")
#' 
#' # Plot the distribution of days between tentatives
#' plot(fML, result="Abort", xlim=c(0, 15))
#' #' 
#' ######### Fit using ML and non parametric ECF
#' 
#' par <- c(ECF.2 = 0.044151921569961131, 
#'          ECF.3 = 2.0020778325280748, 
#'          ECF.4 = 2.6128345101617083, 
#'          ECF.5 = 2.6450582416622375, 
#'          ECF.6 = 2.715198206774927, 
#'          ECF.7 = 2.0288031327239904, 
#'          ECF.8 = 1.0028041546528881, 
#'          ECF.9 = 0.70977432157689235, 
#'          ECF.10 = 0.086052204035003091, 
#'          ECF.11 = 0.011400419961702518, 
#'          ECF.12 = 0.001825219438794076, 
#'          ECF.13 = 0.00029398731859899116, 
#'          ECF.14 = 0.002784886479846703, 
#'          meanIP = 9.9887100433529721, 
#'          sdIP = 0.10580250625108811, 
#'          minIP = 6.5159124624132048, 
#'          pAbort = 2.5702251748938956, 
#'          meanAbort = 2.2721679285648841, 
#'          sdAbort = 0.52006431730489933, 
#'          pCapture = 0.079471782729506113)
#'          
#' fML_NP <- IPFit(x=par, 
#'              fixed.parameters=c(N=20000),
#'              data=data, 
#'              verbose=FALSE, 
#'              model="ML")
#'              
#' par <- fML_NP$ML$par
#' 
#' fML_NP <- IPFit(x=par, 
#'              fixed.parameters=c(N=1000000),
#'              data=data, 
#'              verbose=FALSE, 
#'              model="ML")
#'              
#' par <- c(ECF.2 = 0.016195025683080871, 
#'          ECF.3 = 2.0858089267994315, 
#'          ECF.4 = 3.1307578727979348, 
#'          ECF.5 = 2.7495760827322622, 
#'          ECF.6 = 2.8770821670450939, 
#'          ECF.7 = 2.1592708144943145, 
#'          ECF.8 = 1.0016227335391867, 
#'          ECF.9 = 0.80990178270345259, 
#'          ECF.10 = 0.081051214954249967, 
#'          ECF.11 = 0.039757901443389344, 
#'          ECF.12 = 6.3324056808464527e-05, 
#'          ECF.13 = 0.00037500864146146936, 
#'          ECF.14 = 0.0010383506745475582, 
#'          meanIP = 10.004121090603523, 
#'          sdIP = 0.10229422354470977, 
#'          minIP = 6.5051758088487883, 
#'          pAbort = 2.5335985958484839, 
#'          meanAbort = 2.3145895392189173, 
#'          sdAbort = 0.51192514362374153, 
#'          pCapture = 0.055440514236842105, 
#'          DeltameanIP = -0.046478049165483697)
#' 
#' fML_NP_Delta <- IPFit(x=par, 
#'              fixed.parameters=c(N=20000),
#'              data=data, 
#'              verbose=FALSE, 
#'              model="ML")
#'              
#' par <- fML_NP_Delta$ML$par
#'              
#' fML_NP_Delta <- IPFit(x=par, 
#'              fixed.parameters=c(N=1000000),
#'              data=data, 
#'              verbose=FALSE, 
#'              model="ML")
#'              
#' # Test for stability of -Ln L value according to N
#' grandL.mean <- NULL
#' grandL.sd <- NULL
#' N <- c(10000, 20000, 30000, 40000, 50000,
#'             60000, 70000, 80000, 90000,  
#'             100000, 200000, 300000, 400000, 500000, 
#'             600000, 700000, 800000, 900000,  
#'             1000000)
#' for (Ni in N) {
#'     print(Ni)
#'     smallL <- NULL
#'     for (replicate in 1:100) {
#'          smallL <- c(smallL, 
#'          getFromNamespace(".IPlnL", ns="phenology")
#'                  (x=par, fixed.parameters=c(N=Ni), data=data))
#'     }
#'     grandL.mean <- c(grandL.mean, mean(smallL))
#'     grandL.sd <- c(grandL.sd, sd(smallL))
#' }
#' 
#' grandL.mean <- c(242.619750064524, 239.596145944548, 238.640010536147, 237.965573853263, 
#' 237.727506424543, 237.240740566494, 237.527948232993, 237.297225856515, 
#' 237.17073080938, 237.103397800143, 236.855939567838, 
#' 236.704861853456, 236.82264801458, 236.606065021519, 236.685930841831, 
#' 236.697562908131, 236.568003663293, 236.58097471402, 236.594282543024
#' )
#' grandL.sd <- c(6.54334049298099, 3.04916614991682, 2.57932397492509, 2.15990307710982, 
#' 1.59826856034413, 1.54505295915354, 1.59734964880484, 1.41845032728396, 
#' 1.43096821211286, 1.20048923027244, 0.912467350448495, 
#' 0.75814052890774, 0.668841336554019, 0.539505594152166, 0.554662419326559, 
#' 0.501551009304687, 0.415199780254872, 0.472274287714195, 0.386237047201706
#' )
#' 
#' plot_errbar(x=N, y=grandL.mean, errbar.y = 2*grandL.sd, 
#'             xlab="N", ylab="-Ln L (2 SD)", bty="n", las=1)
#'              
#' # Plot the fitted ECF
#' plot(fML_NP_Delta, result="ECF")
#' 
#' # Plot the Internesting Period distribution
#' plot(fML_NP_Delta, result="IP")
#' 
#' # Plot the distribution of days between tentatives
#' plot(fML_NP_Delta, result="Abort", xlim=c(0, 15))
#' 
#' print(paste("Probability of capture", invlogit(-fML_NP_Delta$ML$par["pCapture"])))
#' # Confidence interval at 95%
#' print(paste(invlogit(-fML_NP_Delta$ML$par["pCapture"]-1.96*fML_NP_Delta$ML$SE["pCapture"]), "-", 
#' invlogit(-fML_NP_Delta$ML$par["pCapture"]+1.96*fML_NP_Delta$ML$SE["pCapture"])))
#' 
#' print(paste("Probability of abort", invlogit(-fML_NP_Delta$ML$par["pAbort"])))
#' # Confidence interval at 95%
#' print(paste(invlogit(-fML_NP_Delta$ML$par["pAbort"]-1.96*fML_NP_Delta$ML$SE["pAbort"]), "-", 
#' invlogit(-fML_NP_Delta$ML$par["pAbort"]+1.96*fML_NP_Delta$ML$SE["pAbort"])))
#'              
#' compare_AIC(parametric=fML$ML, 
#'             nonparameteric=fML_NP$ML, 
#'             nonparametricDelta=fML_NP_Delta$ML)
#' 
#' ######### Fit using Metropolis-Hastings algorithm
#' # ECF.1 = 1 is fixed
#' par <- c(ECF.2 = 0.044151921569961131, 
#'          ECF.3 = 2.0020778325280748, 
#'          ECF.4 = 2.6128345101617083, 
#'          ECF.5 = 2.6450582416622375, 
#'          ECF.6 = 2.715198206774927, 
#'          ECF.7 = 2.0288031327239904, 
#'          ECF.8 = 1.0028041546528881, 
#'          ECF.9 = 0.70977432157689235, 
#'          ECF.10 = 0.086052204035003091, 
#'          ECF.11 = 0.011400419961702518, 
#'          ECF.12 = 0.001825219438794076, 
#'          ECF.13 = 0.00029398731859899116, 
#'          ECF.14 = 0.002784886479846703, 
#'          meanIP = 9.9887100433529721, 
#'          sdIP = 0.10580250625108811, 
#'          minIP = 6.5159124624132048, 
#'          pAbort = 2.5702251748938956, 
#'          meanAbort = 2.2721679285648841, 
#'          sdAbort = 0.52006431730489933, 
#'          pCapture = 0.079471782729506113)
#' 
#' df <- data.frame(Density=rep("dunif", length(par)), 
#' Prior1=c(rep(0, 13), 8, 0.001, 0, -8, 0, 0.001, -8), 
#' Prior2=c(rep(10, 13), 12, 1, 10, 8, 2, 1, 8), 
#' SDProp=unname(c(ECF.2 = 6.366805760909012e-05, 
#'                 ECF.3 = 6.366805760909012e-05, 
#'                 ECF.4 = 6.366805760909012e-05, 
#'                 ECF.5 = 6.366805760909012e-05, 
#'                 ECF.6 = 6.366805760909012e-05, 
#'                 ECF.7 = 6.366805760909012e-05, 
#'                 ECF.8 = 6.366805760909012e-05, 
#'                 ECF.9 = 6.366805760909012e-05, 
#'                 ECF.10 = 6.366805760909012e-05, 
#'                 ECF.11 = 6.366805760909012e-05, 
#'                 ECF.12 = 6.366805760909012e-05, 
#'                 ECF.13 = 6.366805760909012e-05, 
#'                 ECF.14 = 6.366805760909012e-05, 
#'                 meanIP = 6.366805760909012e-05, 
#'                 sdIP = 6.366805760909012e-05, 
#'                 minIP = 6.366805760909012e-05, 
#'                 pAbort = 6.366805760909012e-05, 
#'                 meanAbort = 6.366805760909012e-05, 
#'                 sdAbort = 6.366805760909012e-05, 
#'                 pCapture = 6.366805760909012e-05)),               
#' Min=c(rep(0, 13), 8, 0.001, 0, -8, 0, 0.001, -8), 
#' Max=c(rep(10, 13), 12, 1, 10, 8, 2, 1, 8),
#' Init=par, stringsAsFactors = FALSE)
#' rownames(df)<- names(par)
#' 
#' fMH <- IPFit(parametersMH=df, 
#' fixed.parameters=c(N=10000),
#' data=data, 
#' verbose=FALSE, 
#' n.iter = 10000,
#' n.chains = 1, n.adapt = 100, thin = 1, trace = TRUE,
#' adaptive = TRUE, 
#' model="MH")
#' 
#' # Plot the fitted ECF
#' plot(fMH, result="ECF")
#' 
#' }
#' @export

IPFit <- function (x = NULL, fixed.parameters = NULL, 
                   data = stop("Formated data must be provided"), 
          method = c("Nelder-Mead", "BFGS"), 
          control = list(trace = 1,  REPORT = 100, maxit = 500), 
          itnmax = c(500, 100), hessian = TRUE, 
          verbose = TRUE, parallel = TRUE, model = c("MH", "ML"), parametersMH, 
          n.iter = 10000, n.chains = 1, n.adapt = 100, thin = 30, trace = TRUE, 
          adaptive = TRUE, adaptive.lag = 500, 
          adaptive.fun = function(x) {
            ifelse(x > 0.234, 1.3, 0.7)
          }, intermediate = NULL, filename = "intermediate.Rdata") 
{
  resultML <- list()
  resultMH <- list()
  IPlnL <- getFromNamespace(".IPlnL", ns = "phenology")
  if (any(model == "ML")) {
    if (any(names(x) == "ECF.1")) {
      warning("The parameter ECF.1 is removed from the list of parameter. It must be fixed to 1.")
      x <- x[-which(names(x) == "ECF.1")]
    }
    repeat {
      o <- try(suppressWarnings(optimx::optimx(par = x, 
                                               data = data, fixed.parameters = fixed.parameters, 
                                               fn = IPlnL, method = method, itnmax = itnmax, 
                                               control = modifyList(control, list(dowarn = FALSE, 
                                                                                  follow.on = TRUE, kkt = FALSE)), hessian = FALSE, 
                                               parallel = parallel, verbose = verbose)), silent = TRUE)
      minL <- nrow(o)
      nm <- names(x)
      colnames(o)[1:length(nm)] <- nm
      x <- unlist(o[minL, nm])
      conv <- o[minL, "convcode"]
      value <- o[minL, "value"]
      if (any(names(x) == "meanECF")) 
        x["meanECF"] <- abs(x["meanECF"])
      if (any(substr(names(x), 1, 3) == "ECF")) 
        x[substr(names(x), 1, 3) == "ECF"] <- abs(x[substr(names(x), 
                                                           1, 3) == "ECF"])
      if (any(names(x) == "sdCF")) 
        x["sdCF"] <- abs(x["sdCF"])
      if (any(names(x) == "meanIP")) 
        x["meanIP"] <- abs(x["meanIP"])
      if (any(names(x) == "sdIP")) 
        x["sdIP"] <- abs(x["sdIP"])
      if (any(names(x) == "minIP")) 
        x["minIP"] <- abs(x["minIP"])
      if (any(names(x) == "meanAbort")) 
        x["meanAbort"] <- abs(x["meanAbort"])
      if (any(names(x) == "sdAbort")) 
        x["sdAbort"] <- abs(x["sdAbort"])
      if (conv == 0) 
        break
      message("Convergence is not achieved. Optimization continues !")
    }
    resultML$par <- x
    resultML$value <- value
    resultML$convergence <- conv
    resultML$fixed.parameters <- fixed.parameters
    if (hessian) {
      if (!requireNamespace("numDeriv", quietly = TRUE)) {
        stop("numDeriv package is absent; Please install it first")
      }
      message("Estimation of the standard error of parameters. Be patient please.")
      mathessian <- try(getFromNamespace("hessian", ns = "numDeriv")(func = IPlnL, 
                                                                     fixed.parameters = fixed.parameters, data = data, 
                                                                     x = x, method = "Richardson"), silent = TRUE)
      if (substr(mathessian[1], 1, 5) == "Error") {
        res_se <- rep(NA, length(x))
        names(res_se) <- names(x)
      }      else {
        rownames(mathessian) <- colnames(mathessian) <- names(x)
        resultML$hessian <- mathessian
        res_se <- SEfromHessian(mathessian)
      }
    }    else {
      warning("Standard errors are not estimated.")
      mathessian <- NULL
      res_se <- rep(NA, length(x))
      names(res_se) <- names(x)
    }
    resultML$SE <- res_se
    resultML$AIC <- 2 * resultML$value + 2 * length(x)
    resultML$data <- data
    par <- x
    pfixed <- fixed.parameters
  }
  if (any(model == "MH")) {
    resultMH <- MHalgoGen(likelihood = IPlnL, parameters = parametersMH, 
                          fixed.parameters = fixed.parameters, data = data, 
                          verbose = verbose, parallel = parallel, n.iter = n.iter, 
                          n.chains = n.chains, n.adapt = n.adapt, thin = thin, 
                          trace = trace, adaptive = adaptive, adaptive.lag = adaptive.lag, 
                          adaptive.fun = adaptive.fun, intermediate = intermediate, 
                          filename = filename)
    par <- as.parameters(resultMH)
    pfixed <- fixed.parameters
  }
  model <- IPModel(par = c(par, pfixed))
  result <- list(ML = resultML, MH = resultMH, model = model)
  result <- addS3Class(result, "IP")
  return(result)
}
