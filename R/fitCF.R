#' fitCF fit a model of Clutch Frequency for marine turtles.
#' @title Fit a model of Clutch Frequency for marine turtles.
#' @author Marc Girondot
#' @return Return a list of class ECFOCF with the fit information.\cr
#' The list has the following items:\cr
#' \itemize{
#'   \item \code{data}: The observations to be fitted
#'   \item \code{par}: The fitted parameters
#'   \item \code{SE}: The standard error of parameters if hessian is TRUE
#'   \item \code{value}: The -log likelihood of observations within the fitted model
#'   \item \code{AIC}: The AIC of fitted model
#'   \item \code{mu}: The vector of fitted mu values
#'   \item \code{sd}: The vector of fitted sd values
#'   \item \code{prob}: The vector of fitted capture probabilities
#'   \item \code{a}: The vector of fitted capture probabilities multiplier
#'   \item \code{OTN}: The vector of fitted relative probabilities of contribution
#'   \item \code{period_categories}: A list with the different period probabilities as named vectors for each category
#'   \item \code{period}: The combined period probabilities using OTN as named vector
#'   \item \code{CF_categories}: A list with the different CF probabilities as named vectors for each category
#'   \item \code{CF}: The combined CF probabilities  using OTN as named vector
#'   \item \code{ECFOCF_categories}: A list with the different probability ECFOCF tables for each category
#'   \item \code{ECFOCF}: The combined table of ECFOCF  using OTN probabilities tables
#'   \item \code{ECFOCF_0}: The combined table of ECFOCF probabilities tables  using OTN without the OCF=0
#'   \item \code{SE_df}: A data.frame with SE and 95\% confidence intervals for meanx and vx (mean and variance of clutch frequency for x category), OTNx (proportion for x category), and probx (capture probability for x category)
#' }
#' @param x Initial parameters to be fitted
#' @param fixed.parameters Parameters that are fixed.
#' @param data CMR data formated using TableECFOCF()
#' @param method Method to be used by optimx()
#' @param itnmax A vector with maximum iterations for each method.
#' @param control List of controls for optimx()
#' @param hessian Logical to estimate SE of parameters
#' @param parallel If TRUE, will use parallel computing for ECFOCF_f()
#' @param verbose If TRUE, print the parameters at each step
#' @description This function fits a model of clutch frequency.\cr
#' This model is an enhanced version of the one published by Briane et al. (2007).\cr
#' Parameters are \code{mu} and \code{sd} being the parameters of a  
#' distribution used to model the clutch frequency.\cr
#' This distribution is used only as a guide but has not statistical meaning.\cr
#' The parameter \code{p} is the -logit probability that a female is seen 
#' on the beach for a particular nesting event. It includes both the probability 
#' that it is captured but also the probability that it uses that specific beach.\cr
#' Several categories of females can be included in the model using index after 
#' the name of the parameter, for example \code{mu1}, \code{sd1} and \code{mu2}, 
#' \code{sd2} indicates that two categories of females with different clutch 
#' frequencies distribution are present. Similarly \code{p1} and \code{p2} indicates 
#' that two categories of females with different capture probabilities are present.\cr
#' If more than one category is used, then it is necessary to include the 
#' parameter \code{OTN} to indicate the relative frequencies of each category. 
#' If two categories are used, one \code{OTN} parameter named \code{ONT1} must 
#' be included. The \code{OTN2} is forced to be 1. Then the relative frequency 
#' for category 1 is \code{OTN1/(OTN1+1)} and for category 2 is \code{1/(OTN1+1)}. 
#' Same logic must be applied for 3 and more categories with always the last one 
#' being fixed to 1.\cr
#' if p or a are equal to -Inf, the probability is 0 and if they are equal to 
#' +Inf, the probability is 1.\cr
#' The best way to indicate capture probability for 3D model (OCF, ECF, Period) 
#' is to indicate p.period common for all categories and a1, a2, etc for each category. 
#' The capture probability for category 1 will be p.period * a1, and for category 2 
#' will be p.period * a2, etc. 
#' In this case, the parameters p.period should be indicated in fitted parameters 
#' as well as a1, but a2 must be fixed to +Inf in fixed.parameters. Then the capture 
#' probability for category 2 will be p.period and for category 1 a1 * p.period.
#' @family Model of Clutch Frequency
#' @seealso Briane J-P, Rivalan P, Girondot M (2007) The inverse problem applied 
#'             to the Observed Clutch Frequency of Leatherbacks from Yalimapo beach, 
#'             French Guiana. Chelonian Conservation and Biology 6:63-69
#' @seealso Fossette S, Kelle L, Girondot M, Goverse E, Hilterman ML, Verhage B, 
#'          Thoisy B, de, Georges J-Y (2008) The world's largest leatherback 
#'          rookeries: A review of conservation-oriented research in French 
#'          Guiana/Suriname and Gabon. Journal of Experimental Marine Biology 
#'          and Ecology 356:69-82
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
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
#' # Non parametric model for clutch frequency
#' o_mu1p1_CFnp <- fitCF(x = c(mu.1 = 18.246619595610383, 
#'                        mu.2 = 4.2702163522832892, 
#'                        mu.3 = 2.6289986859556458, 
#'                        mu.4 = 3.2496360919228611, 
#'                        mu.5 = 2.1602522716550943, 
#'                        mu.6 = 0.68617023351032846, 
#'                        mu.7 = 4.2623607001877026, 
#'                        mu.8 = 1.1805600042630455, 
#'                        mu.9 = 2.2786176350939731, 
#'                        mu.10 = 0.47676265496204945, 
#'                        mu.11 = 5.8988238539197062e-08, 
#'                        mu.12 = 1.4003187851424953e-07, 
#'                        mu.13 = 2.4128444894899776e-07, 
#'                        mu.14 = 2.4223748020049825e-07, 
#'                        p = 0.32094401970037578), 
#'                  fixed.parameters=c(mu.15 = 1E-10), 
#'                  data=ECFOCF_2002, hessian = TRUE)
#'                  
#' o_mu2p1 <- fitCF(x = c(mu1 = 1.2190766766978423, 
#'                      sd1 = 0.80646454821956925, 
#'                      mu2 = 7.1886819592223246, 
#'                      sd2 = 0.18152887523015518, 
#'                      p = 0.29347220802963259, 
#'                      OTN = 2.9137627675219533), 
#'                   fixed.parameters=NULL,
#'                   data=ECFOCF_2002, hessian = TRUE)
#' 
#' o_mu1p2 <- fitCF(x = c(mu = 5.3628701816871462, 
#'                      sd = 0.39390555498088764, 
#'                      p1 = 0.61159637544418755, 
#'                      p2 = -2.4212753004659189, 
#'                      OTN = 0.31898004668901009),
#'                  data=ECFOCF_2002, hessian = TRUE)
#'                  
#' o_mu2p2 <- fitCF(x = c(mu1 = 0.043692606004492131, 
#'                    sd1 = 1.9446036983033428, 
#'                    mu2 = 7.3007868915644751, 
#'                    sd2 = 0.16109296152913491, 
#'                    p1 = 1.6860260469536992, 
#'                    p2 = -0.096816113083788985, 
#'                    OTN = 2.2604431232973501), 
#'                   data=ECFOCF_2002, hessian = TRUE)
#'
#' compare_AIC(mu1p1=o_mu1p1_CFp, 
#'             mu2p1=o_mu2p1, 
#'             mu1p2=o_mu1p2, 
#'             mu2p2=o_mu2p2)
#'                  
#' o_mu3p3 <- fitCF(x = c(mu1 = 0.24286312214288761, 
#'                             sd1 = 0.34542255091729313, 
#'                             mu2 = 5.0817174343025551, 
#'                             sd2 = 1.87435099405695, 
#'                             mu3 = 5.2009265101740683, 
#'                             sd3 = 1.79700447678357, 
#'                             p1 = 8.8961708614726156, 
#'                             p2 = 0.94790116453886453, 
#'                             p3 = -0.76572930634505421, 
#'                             OTN1 = 1.2936848663276974, 
#'                             OTN2 = 0.81164278235645926),
#'                  data=ECFOCF_2002, hessian = TRUE)
#'                  
#' 
#' o_mu3p1 <- fitCF(x = structure(c(0.24387978183477, 
#'                                    1.2639261745506, 
#'                                    4.94288464711349, 
#'                                    1.945082889758, 
#'                                    4.9431672350811, 
#'                                    1.287663104591, 
#'                                    0.323636536050397, 
#'                                    1.37072039291397, 
#'                                    9.28055412564559e-06), 
#'                                   .Names = c("mu1", "sd1", "mu2", 
#'                                              "sd2", "mu3", "sd3", 
#'                                              "p", "OTN1", "OTN2")),
#'                  data=ECFOCF_2002, hessian = TRUE)
#'                  
#' 
#' o_mu1p3 <- fitCF(x = structure(c(4.65792402108387, 
#'                                    1.58445909785, 
#'                                    -2.35414198317177, 
#'                                    0.623757854800649, 
#'                                    -3.62623634029326, 
#'                                    11.6950204755787, 
#'                                    4.05273728846523), 
#'                                    .Names = c("mu", "sd", 
#'                                               "p1", "p2", "p3", 
#'                                               "OTN1", "OTN2")),
#'                  data=ECFOCF_2002, hessian = TRUE)
#'                  
#' compare_AIC(mu1p1=o_mu1p1, 
#'             mu2p1=o_mu2p1, 
#'             mu1p2=o_mu1p2, 
#'             mu2p2=o_mu2p2, 
#'             mu3p3=o_mu3p3, 
#'             mu1p3=o_mu1p3, 
#'             mu3p1=o_mu3p1)
#'             
#'  # 3D model for (ECF, OCF, period)           
#'             
#' ECFOCF_2002 <- TableECFOCF(MarineTurtles_2002, 
#'                            date0=as.Date("2002-01-01"))
#' 
#' fp <- rep(0, dim(ECFOCF_2002)[3])
#' names(fp) <- paste0("p.", formatC(1:(dim(ECFOCF_2002)[3]), width=2, flag="0"))
#' par <- c(mu = 2.6404831115214353, 
#'         sd = 0.69362774786433479, 
#'         mu_season = 12.6404831115214353, 
#'         sd_season = 1.69362774786433479)
#' par <- c(par, fp[attributes(ECFOCF_2002)$table["begin"]:
#'                  attributes(ECFOCF_2002)$table["end"]])
#' # The value of p (logit -capture probability) out of the period 
#' # of monitoring is set to +Inf (capture probability=1)
#' # to indicate that no turtle is nesting in the period out of 
#' # monitoring time
#' # p is set to -Inf (capture probability=0) to indicate that no
#' # monitoring has been done but some turtles could have been present.
#' fixed.parameters <- c(p=+Inf)
#' # The fitted values are:
#' par <- c(mu = 2.4911638591178051, 
#'          sd = 0.96855483039640977, 
#'          mu_season = 13.836059118657793, 
#'          sd_season = 0.17440085345943984, 
#'          p.10 = 1.3348233607728222, 
#'          p.11 = 1.1960387774393837, 
#'          p.12 = 0.63025680979544774, 
#'          p.13 = 0.38648155002707452, 
#'          p.14 = 0.31547864054366048, 
#'          p.15 = 0.19720001827017075, 
#'          p.16 = 0.083199496372073328, 
#'          p.17 = 0.32969130595897905, 
#'          p.18 = 0.36582777525265819, 
#'          p.19 = 0.30301248314170637, 
#'          p.20 = 0.69993987591518514, 
#'          p.21 = 0.13642423871641118, 
#'          p.22 = -1.3949268190534629)
#' 
#' o_mu1p1season1 <- fitCF(x=par, data=ECFOCF_2002, 
#'                         fixed.parameters=fixed.parameters)
#'
#' # Same model but with two different models of capture probabilities
#'                         
#' fp <- rep(0, dim(ECFOCF_2002)[3])
#' names(fp) <- paste0("p1.", formatC(1:(dim(ECFOCF_2002)[3]), width=2, flag="0"))
#' par <- c(mu = 2.6404831115214353, 
#'         sd = 0.69362774786433479, 
#'         mu_season = 12.6404831115214353, 
#'         sd_season = 1.69362774786433479)
#' par <- c(par, fp[attributes(ECFOCF_2002)$table["begin"]:
#'                  attributes(ECFOCF_2002)$table["end"]])
#' names(fp) <- paste0("p2.", formatC(1:(dim(ECFOCF_2002)[3]), width=2, flag="0"))
#' par <- c(par, fp[attributes(ECFOCF_2002)$table["begin"]:
#'                  attributes(ECFOCF_2002)$table["end"]])
#' fixed.parameters <- c(p1=+Inf, p2=+Inf)
#' 
#' o_mu1p2season1 <- fitCF(x=par, data=ECFOCF_2002, 
#'                         fixed.parameters=fixed.parameters)
#' 
#' # Here the two different capture probabilities are different 
#' # by a constant:
#' # p1=invlogit(-p)     [Note that invlogit(-a1) = 1]
#' # p2=invlogit(-p)*invlogit(-a2)
#' 
#' fp <- rep(0, dim(ECFOCF_2002)[3])
#' names(fp) <- paste0("p.", formatC(1:(dim(ECFOCF_2002)[3]), width=2, flag="0"))
#' par <- c(mu = 2.6404831115214353, 
#'         sd = 0.69362774786433479, 
#'         mu_season = 12.6404831115214353, 
#'         sd_season = 1.69362774786433479, 
#'         a2=0)
#' par <- c(par, fp[attributes(ECFOCF_2002)$table["begin"]:
#'                  attributes(ECFOCF_2002)$table["end"]])
#' fixed.parameters <- c(a1=+Inf, p=+Inf)
#' 
#' o_mu1p1aseason1 <- fitCF(x=par, data=ECFOCF_2002, 
#'                         fixed.parameters=fixed.parameters)
#'                         data=ECFOCF_2002)
#'  
#' }
#' @export

# Lancement du fit ####

# library("phenology");load(file="/Users/marcgirondot/Documents/Espace_de_travail_R/Remigration/CF_R/dataOut/fit2002_CF.Rdata"); for (i in names(fit2002_CF)) assign(paste0("o_", i), fit2002_CF[[i]])

fitCF <- function(x=c(mu=4, sd=100, p=0),
                  fixed.parameters=NULL, 
                  data=stop("Data formated with TableECFOCF() must be provided"),
                  method = c("Nelder-Mead","BFGS"), 
                  control=list(trace=1, REPORT=100, maxit=500),
                  itnmax=c(500, 100), 
                  hessian=TRUE, parallel=TRUE, verbose=FALSE) {
  
  
#  x=c(mu=4, sd=100, p=-1);
#  fixed.parameters=NULL;
#  data=NULL;
#  method = c("Nelder-Mead","BFGS");
#  control=list(trace=1, REPORT=100, maxit=500);
#  itnmax=c(500, 100);
#  hessian=TRUE; parallel=TRUE
  
  MaxNests <- max(dim(data)[c(1, 2)])-1
  
  repeat {
    o <- try(suppressWarnings(optimx::optimx(par = x,
                            data=data, 
                            fixed.parameters=fixed.parameters,
                            fn=lnLCF, 
                            method=method, 
                            itnmax=itnmax, 
                            control=modifyList(control, list(dowarn=FALSE, follow.on=TRUE, kkt=FALSE)), 
                            hessian=FALSE, parallel=parallel, verbose=verbose)), silent=TRUE)
    
    minL <- nrow(o)
    nm <- names(x)
    colnames(o)[1:length(nm)] <- nm
    # nm <- gsub("-", ".", nm)
    # 
    x <- unlist(o[minL, nm])
    conv <- o[minL, "convcode"]
    value <- o[minL, "value"]
    
    x[substr(names(x), 1, 2)=="mu"] <- abs(x[substr(names(x), 1, 2)=="mu"])
    x[substr(names(x), 1, 2)=="sd"] <- abs(x[substr(names(x), 1, 2)=="sd"])
    x[substr(names(x), 1, 3)=="OTN"] <- abs(x[substr(names(x), 1, 3)=="OTN"])
    # C'est déjà inclut dans le mu et le sd
    # x[substr(names(x), 1, 9)=="mu_season"] <- abs(x[substr(names(x), 1, 9)=="mu_season"])
    # x[substr(names(x), 1, 9)=="sd_season"] <- abs(x[substr(names(x), 1, 9)=="sd_season"])
    
    if (conv == 0) break
    # par <- x
    message("Convergence is not achieved. Optimization continues !")
  }
  
  result <- list()
  result$par <- x
  result$value <- value
  result$convergence <- conv
  result$fixed.parameters <- fixed.parameters
  
  if (hessian) {
    
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      stop("numDeriv package is absent; Please install it first")
    }
    
    message("Estimation of the standard error of parameters. Be patient please.")
    
    mathessian <- try(getFromNamespace("hessian", ns="numDeriv")(func=lnLCF, 
                                                                 fixed.parameters=fixed.parameters, 
                                                                 data=data, 
                                                                 x=x, 
                                                                 method="Richardson"
    )
    , silent=TRUE)
    if (substr(mathessian[1], 1, 5)=="Error") {
      res_se <- rep(NA, length(x))
      names(res_se) <- names(x)
    } else {
      
      rownames(mathessian) <- colnames(mathessian) <- names(x)
      result$hessian <- mathessian
      res_se <- SEfromHessian(mathessian)
    }
  } else {
    warning("Standard errors are not estimated.")
    mathessian <- NULL
    res_se <- rep(NA, length(x))
    names(res_se) <- names(x)
  }
  
  result$SE <- res_se
  
  
  result$AIC <- 2*result$value+2*length(x)
  result$data <- data
  
  totx <- c(x, fixed.parameters)
  
  ml <- floor(as.numeric(gsub("[a-zA-z]+", "", names(totx ))))
  if ((length(ml) == 1) | all(is.na(ml)) | (max(c(0, ml), na.rm=TRUE)==0)) {
    mln <- 1
  } else {
    mln <- max(ml, na.rm=TRUE)
  }
  
  
  p <- totx[substr(names(totx), 1, 1)=="p"]
  if (length(p)>1) {
    np <- gsub("p([0-9\\.]*)", "\\1", names(p))
    np[np == ""] <- "0"
    p <- p[order(as.numeric(np))]
  }
  
  a <- totx[substr(names(totx), 1, 1)=="a"]
  if (identical(a, structure(numeric(0), .Names = character(0)))) {
    a <- rep(Inf, mln) # Vaudra 1
    names(a) <- paste0("a", as.character(1:mln))
  }
  if (length(a)>1) a <- a[order(as.numeric(gsub("a([0-9]+)", "\\1", names(a))))]
  
  mu <- totx[(substr(names(totx), 1, 2)=="mu") & (substr(names(totx), 1, 9)!="mu_season")]
  if (length(mu)>1) mu <- mu[order(as.numeric(gsub("mu([0-9\\.]+)", "\\1", names(mu))))]
  sd <- totx[(substr(names(totx), 1, 2)=="sd") & (substr(names(totx), 1, 9)!="sd_season")]
  if (identical(sd, structure(numeric(0), .Names = character(0)))) sd <- c(sd=NA)
  if (length(sd)>1) sd <- sd[order(as.numeric(gsub("sd([0-9]+)", "\\1", names(sd))))]
  mu_season <- totx[substr(names(totx), 1, 9)=="mu_season"]
  if (length(mu_season)>1) mu_season <- mu_season[order(as.numeric(gsub("mu_season([0-9]+)", "\\1", names(mu_season))))]
  
  sd_season <- totx[substr(names(totx), 1, 9)=="sd_season"]
  if (length(sd_season)>1) sd_season <- sd_season[order(as.numeric(gsub("sd_season([0-9]+)", "\\1", names(sd_season))))]
  
  if (mln>1) {
    OTN <- abs(totx[substr(names(totx), 1, 3)=="OTN"])
    if (length(OTN)>1) OTN <- OTN[order(as.numeric(gsub("OTN([0-9]+)", "\\1", names(OTN))))]
    OTN <- c(OTN, 1)
    OTN <- c(OTN, rep(OTN[length(OTN)], mln-length(OTN)))
    names(OTN) <- paste0("OTN", 1:mln)
    OTN <- OTN/sum(OTN)
  } else {
    OTN <- c(OTN1=1)
  }
  
  result$OTN <- OTN
  
  p <- 1/(1+exp(-p))
  a <- 1/(1+exp(-a))
  
  result$prob <- p
  result$a <- a
  
  # mu <- abs(c(mu, rep(mu[length(mu)], mln-length(mu))))
  # names(mu) <- paste0("mu", 1:mln)
  
  if (mln > 1) {
    if (any(names(mu)=="mu")) {
      mu_ref <- mu[names(mu)=="mu"]
      mu_ec <- NULL
      for (i in 1:mln) {
        if (all(!grepl(paste0("mu", i), names(mu)))) {
          mu_ec <- c(mu_ec, structure(unname(mu_ref), .Names=paste0("mu", i)))
        } else {
          mu_ec <- c(mu_ec, mu[grepl(paste0("mu", i), names(mu))])
        }
      }
      mu <- mu_ec
    }
    if (any(grepl("mu\\.+", names(mu)))) {
      mu_ref <- mu[grepl("mu\\.+", names(mu))]
      
      mu_ec <- NULL
      for (i in 1:mln) {
        if (all(!grepl(paste0("mu", i), names(mu)))) {
          mu_ec <- c(mu_ec, structure(unname(mu_ref), .Names=gsub("mu(\\.[0-9]+)", paste0("mu", i, "\\1"), names(mu_ref))))
        } else {
          mu_ec <- c(mu_ec, mu[grepl(paste0("mu", i), names(mu))])
        }
      }
      mu <- mu_ec
    }
  } else {
    names(mu) <- gsub("mu[0-9]*(\\.*[0-9]*)", "mu1\\1", names(mu))
  }
  
  result$mu <- mu
  
  sd <- abs(c(sd, rep(sd[length(sd)], mln-length(sd))))
  names(sd) <- paste0("sd", 1:mln)
  
  result$sd <- sd
  
  if (!identical(unname(mu_season), numeric(0))) {
  
  mu_season <- abs(c(mu_season, rep(mu_season[length(mu_season)], mln-length(mu_season))))
  names(mu_season) <- paste0("mu_season", 1:mln)
  
  result$mu_season <- mu_season
  
  sd_season <- abs(c(sd_season, rep(sd_season[length(sd_season)], mln-length(sd_season))))
  names(sd_season) <- paste0("sd_season", 1:mln)
  
  result$sd_season <- sd_season
  }
  
  OCFECF <- data
  OCFECF[] <- 0
  
  CF <- rep(0, MaxNests)
  names(CF) <- paste0("CF", as.character(1:MaxNests))
  
  if (dim(data)[3] != 1) {
    period <- structure(rep(0, dim(data)[3]-MaxNests), 
                        .Names=paste0("period", formatC(1:(dim(data)[3]-MaxNests), width=2, flag="0")))
  } else {
    period <- NA
  }
  
  OCFECF_categories <- list()
  OCFECF_0_categories <- list()
  CF_categories <- list()
  if (!is.na(period[1])) period_categories <- list()
  
  for (i in 1:mln) {
    nm <- paste0("p", as.character(i))
    
    pvrai <- p[(names(p) == "p") | (substr(names(p), 1, nchar(nm))==nm) | 
                 (substr(names(p), 1, 2)=="p.")]
    
    # Il faut rajouter a[i] seulement si avec une période
    pvrai[grepl("\\.", names(pvrai))] <- pvrai[grepl("\\.", names(pvrai))] * a[paste0("a", i)]
    
    names(pvrai)[substr(names(pvrai), 1, 2) =="p."] <- paste0("p", i, ".", substr(names(pvrai[substr(names(pvrai), 1, 2) =="p."]), 3, 20))
    names(pvrai)[names(pvrai) =="p"] <- paste0("p", i)
    
    OCFECF_int <- ECFOCF_f(mu=mu[grepl(paste0("mu", i), names(mu))],
                           sd=sd[paste0("sd", i)], 
                           p=pvrai, 
                           MaxNests=MaxNests, 
                           mu_season=mu_season[paste0("mu_season", i)], 
                           sd_season=sd_season[paste0("sd_season", i)], 
                           length_season = dim(data)[3]-MaxNests, 
                           parallel=parallel
    ) 
    class(OCFECF_int) <- "TableECFOCF"
    OCFECF_categories <- c(OCFECF_categories, list(OCFECF_int))
    OCFECF <- OCFECF+ OCFECF_int* OTN[paste0("OTN", i)]
    
    OCFECF_int <- OCFECF_int/(1-sum(OCFECF_int[1, 1, ]))
    OCFECF_int[1, 1, ] <- 0
    
    OCFECF_0_categories <- c(OCFECF_0_categories, list(OCFECF_int))
    
    if (!is.na(sd[paste0("sd", i)])) {
      # Ancienne formule
      CF_int <- dlnorm(1:MaxNests, meanlog=log(abs(mu[paste0("mu", i)])), 
                  sdlog=abs(sd[paste0("sd", i)]))
    } else {
      # Nouvelle formule
      # Je dois sortir les mu classés par ordre croissant
      CF_int <- abs(mu[order(as.numeric(gsub("mu[0-9]*\\.", "", names(mu))))])
      if (length(CF_int) < MaxNests) CF_int <- c(CF_int, rep(1E-10, MaxNests - length(CF_int)))
    }
    
    CF_int <- structure(c(CF_int / sum(CF_int)), .Names=paste0("CF", as.character(1:MaxNests)))
    CF_categories <- c(CF_categories, 
                 list(CF_int))
    CF <- CF+ CF_int * OTN[paste0("OTN", i)]
    
    if (!is.na(period[1])) {
      time_int <- dlnorm(1:(dim(data)[3]-MaxNests), 
                         meanlog=log(abs(mu_season[paste0("mu_season", i)])), 
                         sdlog=abs(sd_season[paste0("sd_season", i)]))
      time_int <- structure(c(time_int / sum(time_int)), 
                            .Names=paste0("period", formatC(1:(dim(data)[3]-MaxNests), width=2, flag="0")))
      
      period_categories <- c(period_categories, list(time_int))
      
      period <- period + time_int * OTN[paste0("OTN", i)]
    }
  }
  
  result$ECFOCF_categories <- OCFECF_categories
  result$CF_categories <- CF_categories
  class(OCFECF) <- "TableECFOCF"
  result$ECFOCF <- OCFECF
  result$CF <- CF
  if (!is.na(period[1])) {
    result$period_categories <- period_categories
    result$period <- period
    
    result$length_season <- dim(data)[3]-MaxNests
  }
  
  
  OCFECF <- OCFECF/(1-sum(OCFECF[1, 1, ]))
  OCFECF[1, 1, ] <- 0
  
  result$ECFOCF_0 <- OCFECF
  result$ECFOCF_0_categories <- OCFECF_0_categories
  
  result$MaxNests <- MaxNests
  
  result$categories <- mln
  
  if (hessian & is.element('car', installed.packages()[,1])) {
    vcov <- try(solve(mathessian), silent = TRUE)
    if (class(vcov) == "try-error") {
      message("Error in inverse of Hessian matrix, some standard errors cannot be calculated")
    }
      # mu
      SE_df <- data.frame(Estimate=numeric(), 
                          SE=numeric(), 
                          "2.5 %"=numeric(),
                          "97.5 %"=numeric())
      colnames(SE_df) <- c("Estimate", "SE", "2.5 %", "97.5 %")
      # SE_df_0 <- SE_df
      
      par_mu <- names(x)[(substr(names(x), 1, 2) == "mu") & (substr(names(x), 1, 9) != "mu_season")]
      for (i in seq_along(par_mu)) {
        if (class(vcov) == "try-error") {
          SE_df[nrow(SE_df)+1, ] <- c(eval(parse(text=x[par_mu[i]])), NA, NA, NA)
        } else {
          SE_df[nrow(SE_df)+1, ] <- unlist(car::deltaMethod(x, par_mu[i], vcov.=vcov)[1, c(1:4), drop = TRUE])
        }
      }
      colnames(SE_df) <- c("Estimate", "SE", "2.5 %", "97.5 %")
      rownames(SE_df) <- par_mu
      
      rownames(SE_df) <- gsub("mu", "mean", par_mu)
      SE_df[, "Estimate"] <- SE_df[, "Estimate"] +1 
      SE_df[, "2.5 %"] <- SE_df[, "2.5 %"] +1 
      SE_df[, "97.5 %"] <- SE_df[, "97.5 %"] +1 
      
      rn <- rownames(SE_df)
      par_mu_season <- names(x)[(substr(names(x), 1, 9) == "mu_season")]
      for (i in seq_along(par_mu_season)) {
        if (class(vcov) == "try-error") {
          SE_df[nrow(SE_df)+1,] <- c(eval(parse(text=x[par_mu_season[i]])), NA, NA, NA)
        } else {
          SE_df[nrow(SE_df)+1,] <- unlist(car::deltaMethod(x, par_mu_season[i], vcov.=vcov)[1, c(1:4), drop = TRUE])
        }
      }
      rownames(SE_df) <- c(rn, gsub("mu_", "mean_", par_mu_season))
      
      # OTN
      
      par_OTN <- names(x)[substr(names(x), 1, 3) == "OTN"]
      rn <- rownames(SE_df)
      if (! identical(par_OTN, character(0))) {
        
        for (i in c(par_OTN, "1")) {
          if (class(vcov) == "try-error") {
            
            denom <- paste0("/(1 + ",paste(paste0("x['", par_OTN, "']") , collapse = "+"), ")", collapse="")
            num <- ifelse(i=="1", "1", paste0("x['", i, "']"))
            SE_df[nrow(SE_df)+1, ] <- c(eval(parse(text=paste0(num, denom))), NA, NA, NA)
          } else {
            denom <- paste0("/(1 + ", paste(par_OTN , collapse = "+"), ")", collapse="")
            
            SE_df[nrow(SE_df)+1, ] <- unlist(car::deltaMethod(x, 
                                                                          paste0(i, denom)
                                                                          , vcov.=vcov)[1, c(1:4), drop = TRUE])
          }
        }
        rownames(SE_df) <- c(rn, paste0("OTN", as.character(1:(nrow(SE_df)-length(rn)))))
      }
      
      rn <- rownames(SE_df)
      # p
      par_p <- names(x)[substr(names(x), 1, 1) == "p"]
      par_a_tot <- names(totx)[substr(names(totx), 1, 1) == "a"]
      par_a <- names(x)[substr(names(x), 1, 1) == "a"]
      
      if (! identical(par_a_tot, character(0))) {
        for (j in par_a_tot) {
          for (i in par_p) {
            
            categp <- gsub("p", "", gsub("\\.[0-9]+", "", i))
            catega <- gsub("p", "", gsub("\\.[0-9]+", "", j))
            
            if ((categp == "") | (catega == categp)) {
              
              if (any(j == par_a)) {
                if (class(vcov) == "try-error") {
                  SE_df[nrow(SE_df)+1, ] <- c(eval(parse(text=paste0("1/(1+exp(", x[i], ")) * 1/(1+exp(", x[catega], "))"))), NA, NA, NA)
                  rn <- c(rn, paste0("1/(1+exp(", i, ")) * 1/(1+exp(", catega, "))"))
                } else {
                  SE_df[nrow(SE_df)+1, ] <- unlist(car::deltaMethod(x, 
                                                                                paste0("1/(1+exp(", i, ")) * 1/(1+exp(", catega, "))")
                                                                                , vcov.=vcov)[1, c(1:4), drop = TRUE])
                  rn <- c(rn, paste0("1/(1+exp(", i, ")) * 1/(1+exp(", catega, "))"))
                }
                  } else {
                    if (class(vcov) == "try-error") {
                      SE_df[nrow(SE_df)+1, ] <- c(eval(parse(text=paste0("1/(1+exp(", x[i], "))"))), NA, NA, NA)
                      rn <- c(rn, row.names = paste0("1/(1+exp(", i, "))"))
                    } else {
                      
                      SE_df[nrow(SE_df)+1, ] <- unlist(car::deltaMethod(x, 
                                                                                paste0("1/(1+exp(", i, "))")
                                                                                , vcov.=vcov)[1, c(1:4), drop = TRUE])
                      rn <- c(rn, row.names = paste0("1/(1+exp(", i, "))"))
                    }
                      }
            }
          }
        }
      } else {
        for (i in par_p) {
          if (class(vcov) == "try-error") {
            SE_df[nrow(SE_df)+1, ] <- c(eval(parse(text=paste0("1/(1+exp(", x[i], "))"))), NA, NA, NA)
            rn <- c(rn, row.names = paste0("1/(1+exp(", i, "))"))
          } else {
            SE_df[nrow(SE_df)+1, ] <- unlist(car::deltaMethod(x, 
                                                       paste0("1/(1+exp(", i, "))")
                                                       , vcov.=vcov)[1, c(1:4), drop = TRUE])
            rn <- c(rn, row.names = paste0("1/(1+exp(", i, "))"))
          }
            }
      }
      rnp <- rn
      rnp <- gsub("1/\\(1\\+exp\\(", "", rnp)
      # rnp <- gsub("1/\\(1 \\+ exp\\(", "", rnp)
      rnp <- gsub("\\)", "", rnp)
      rnp <- gsub("p", "prob", rnp)
      rownames(SE_df) <- rnp

      or <- gsub("[A-Za-z_]+([0-9\\.]+)$", "\\1", rownames(SE_df))
      or <- gsub("[A-Za-z_]+$", "", or)
      or <- gsub("^[A-Za-z_]+", "", or)
      or <- gsub("([0-9\\.]+) \\* ([0-9\\.]+)", "\\2\\1", or)
      or <- ifelse (or=="", 0, or)
      # NA dans certains cas 20/1/2018
      SE_df <- SE_df[order(as.numeric(or)), ]
      
      result$SE_df <- SE_df
    
  }
  
  class(result) <- "ECFOCF"
  return(result)
}

