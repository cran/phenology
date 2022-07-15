#' Tagloss_model returns the daily rate of tag loss.
#' @title Return the daily rate of tag loss.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Return the daily rate of tag loss if hessian is null or a data.frame with distribution of daily rate of tag loss if hessian is not null.
#' @param t Time for which values of model must be estimated
#' @param x A Tagloss fitted model
#' @param par Parameters
#' @param Hessian Hessian matrix of parameters
#' @param mcmc A mcmc result
#' @param model The model of parameter to be used, can be 1, 2, L1, L2, R1 or R2
#' @param model_before Function to be used before estimation of daily tagloss rate
#' @param model_after Function to be used after estimation of daily tagloss rate
#' @param method Can be NULL, "delta", "SE", "Hessian", "MCMC", or "PseudoHessianFromMCMC"
#' @param replicates Number of replicates to estimate se of output
#' @description This function compute a model of daily tag loss rate for days t 
#' based on a set of parameters, par or a fitted tag loss model in x.\cr
#' Parameters are described in \code{\link{Tagloss_fit}}.
#' @family Model of Tag-loss
#' @examples
#' \dontrun{
#' library(phenology)
#' 
#' # Example
#' t <- 1:1000
#' par <- c(D1=200, D2D1=100, D3D2=200, 
#'          A=-logit(0.02), B=-logit(0.05), C=-logit(0.07))
#' y <- Tagloss_model(t, par, model="1")
#' plot(x=t, y, type="l")
#' par <- c(D1_1=200, D2D1_1=100, D3D2_1=200, 
#'          A_1=-logit(0.02), B_1=-logit(0.05), C_1=-logit(0.07))
#' y <- Tagloss_model(t, par, model="1")
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, model="1")
#' 
#' # Fig1A in Rivalan et al. 2005 (note an error for a0; a0 must be negative)
#' par <- c(a0=-1E5, a1=-2000, a2=0, a3=2*max(t), a4=0.1)
#' y <- Tagloss_model(t, par)
#' plot(x=t, y, type="l")
#' 
#' # Fig1B in Rivalan et al. 2005
#' par <- c(a0=-0.5, a1=-2000, a2=-0.001, a3=0, a4=0.1)
#' y <- Tagloss_model(t, par)
#' plot(x=t, y, type="l")
#' 
#' # Fig1C in Rivalan et al. 2005
#' par <- c(a0=-1, a1=-6, a2=0, a3=0, a4=0)
#' y <- Tagloss_model(t, par)
#' plot(x=t, y, type="l")
#' 
#' # Fig1D in Rivalan et al. 2005
#' par <- c(a0=-1, a1=-6, a2=0, a3=0, a4=0.1)
#' y <- Tagloss_model(t, par)
#' plot(x=t, y, type="l")
#' 
#' # Fig1E in Rivalan et al. 2005
#' par <- c(a0=-0.1, a1=-10, a2=-0.2, a3=60, a4=0.1)
#' y <- Tagloss_model(t, par)
#' plot(x=t, y, type="l")
#' 
#' # Fig1F in Rivalan et al. 2005
#' par <- c(a0=-0.1, a1=-10, a2=0.2, a3=60, a4=0.1)
#' y <- Tagloss_model(t, par)
#' plot(x=t, y, type="l")
#' 
#' # Example with fitted data
#' data_f_21 <- Tagloss_format(outLR, model="21")
#' # Without the N20 the computing is much faster
#' data_f_21_fast <- subset(data_f_21, subset=(is.na(data_f_21$N20)))
#' par <- c('D1_2' = 49.086835072129126, 
#'          'D2D1_2' = 1065.0992647723231, 
#'          'D3D2_2' = 6.15531475922079, 
#'          'A_2' = 5.2179675647973758, 
#'          'B_2' = 8.0045560376751386, 
#'          'C_2' = 8.4082505219581876, 
#'          'D1_1' = 177.23337287498103, 
#'          'D2D1_1' = 615.42690323741033, 
#'          'D3D2_1' = 2829.0806609455867, 
#'          'A_1' = 28.500118091731551, 
#'          'B_1' = 10.175426055942701, 
#'          'C_1' = 6.9616630417169398)
#' o <- Tagloss_fit(data=data_f_21_fast, fitted.parameters=par)
#' t <- 1:10
#' y <- Tagloss_model(t, o$par, model="1")
#' y <- Tagloss_model(t, x=o, method=NULL, model="1")
#' y <- Tagloss_model(t, x=o, method="Hessian", model="1", replicates=1000)
#' }
#' @export

Tagloss_model <- function(t=NULL                                                      , 
                          par=NULL                                                    , 
                          Hessian=NULL                                                , 
                          mcmc = NULL                                                 ,
                          model_before = NULL                                         , 
                          model_after = NULL                                          , 
                          model=stop("You must specify which tag loss rate you want."), 
                          method=NULL                                                 , 
                          replicates=NULL                                             , 
                          x=NULL                                                      ) {
  
  if (is.null(replicates)) replicates <- 0
  if (is.null(method) | (replicates == 0)) method <- "null"
  method <- tolower(method)
  method <- match.arg(method, choices = c("null", 
                                          "delta", 
                                          "se", 
                                          "hessian", 
                                          "mcmc", 
                                          "pseudohessianfrommcmc"))
  
  if (is.null(x) & is.null(par)) {
    stop("Both par and x cannot be null at the same time")
  }
  
  if (is.null(t)) {
    if (is.null(x)) {
      stop("Both t and x cannot be null at the same time")
    } else {
      t <- 1:x$days.maximum
    }
  }
  
  if (!is.null(x)) {
    if (is.null(par)) par <- c(x$par, x$fixed.par)
    if (is.null(Hessian)) Hessian <- x$hessian
    if (is.null(model_before)) model_before <- x$model_before
    if (is.null(model_after)) model_after <- x$model_after
  }
  
  if (!is.null(model)) {
    model <- toupper(model)
    model <- match.arg(model, choices = c("1", "2", "L1", "L2", "R1", "R2"))
  } else {
    stop("Model cannot be NULL.")
  }
  
  if (method=="hessian" & is.null(Hessian)) {
    stop("Hessian method cannot be used with an empty Hessian matrix.")
  }
  
  if ((method=="mcmc" | method=="pseudohessianfrommcmc") & is.null(mcmc)) {
    stop("MCMC methods cannot be used without MCMC parameter.")
  }
  
  days.maximum <- max(t)
  
  if (method=="null") {
    
    if (!is.null(model_before)) eval(parse(text=model_before), envir= environment())
    
    p1 <- NULL
    p2 <- NULL
    pL1 <- NULL
    pR1 <- NULL
    pL2 <- NULL
    pR2 <- NULL
    
    # Casale 2017 Model 1
    if (any(names(par) == "CasaleModelIc_1")) p1 <- rep(par["CasaleModelIc_1"], length(t))
    if (any(names(par) == "CasaleModelIc_2")) p2 <- rep(par["CasaleModelIc_2"], length(t))
    if (any(names(par) == "CasaleModelIc")) p2 <- p1 <- rep(par["CasaleModelIc"], length(t))
    
    # Casale 2017 Model II
    if (any(names(par) == "CasaleModelIIa0_1")) p1 <- (1-par["CasaleModelIIa4_1"])/(1+exp(par["CasaleModelIIa0_1"]*(par["CasaleModelIIa1_1"]-t))) + par["CasaleModelIIa4_1"]
    if (any(names(par) == "CasaleModelIIa0_2")) p2 <- (1-par["CasaleModelIIa4_2"])/(1+exp(par["CasaleModelIIa0_2"]*(par["CasaleModelIIa1_2"]-t))) + par["CasaleModelIIa4_2"]
    if (any(names(par) == "CasaleModelIIa0")) p1 <- p2 <- (1-par["CasaleModelIIa4"])/(1+exp(par["CasaleModelIIa0"]*(par["CasaleModelIIa1"]-t))) + par["CasaleModelIIa4"]
    
    # Casale 2017 Model III
    if (any(names(par) == "CasaleModelIIIa0_1")) p1 <- par["CasaleModelIIIa4_1"]/(1+exp(par["CasaleModelIIIa0_1"]*(par["CasaleModelIIIa1_1"]-t)))
    if (any(names(par) == "CasaleModelIIIa0_2")) p2 <- par["CasaleModelIIIa4_2"]/(1+exp(par["CasaleModelIIIa0_2"]*(par["CasaleModelIIIa1_2"]-t)))
    if (any(names(par) == "CasaleModelIIIa0")) p1 <- p2 <- par["CasaleModelIIIa4"]/(1+exp(par["CasaleModelIIIa0"]*(par["CasaleModelIIIa1"]-t)))
    
    # Casale 2017 Model IV
    if (any(names(par) == "CasaleModelIVa0_1")) {
      mint <- par["CasaleModelIVa4_1"]/(1+exp(par["CasaleModelIVa2_1"]*(par["CasaleModelIVa3_1"]-t)))
      p1 <- (0.01-mint)/(1+exp(par["CasaleModelIVa0_1"]*(par["CasaleModelIVa1_1"]-t))) + mint
    }
    if (any(names(par) == "CasaleModelIVa0_2")) {
      mint <- par["CasaleModelIVa4_2"]/(1+exp(par["CasaleModelIVa2_2"]*(par["CasaleModelIVa3_2"]-t)))
      p2 <- (0.01-mint)/(1+exp(par["CasaleModelIVa0_2"]*(par["CasaleModelIVa1_2"]-t))) + mint
    }
    if (any(names(par) == "CasaleModelIVa0")) {
      mint <- par["CasaleModelIVa4"]/(1+exp(par["CasaleModelIVa2"]*(par["CasaleModelIVa3"]-t)))
      p1 <- p2 <- (0.01-mint)/(1+exp(par["CasaleModelIVa0"]*(par["CasaleModelIVa1"]-t))) + mint
    }    
    # Casale 2017 Model 5
    if (any(names(par) == "CasaleModelVa0_1")) {
      maxt <- (0.01-par["CasaleModelVa4_1"])/(1+exp(par["CasaleModelVa2_1"]*(par["CasaleModelVa3_1"]-t))) + par["CasaleModelVa4_1"]
      p1 <- maxt/(1+exp(par["CasaleModelVa0_1"]*(par["CasaleModelVa1_1"]-t)))
    }
    if (any(names(par) == "CasaleModelVa0_2")) {
      maxt <- (0.01-par["CasaleModelVa4_2"])/(1+exp(par["CasaleModelVa2_2"]*(par["CasaleModelVa3_2"]-t))) + par["CasaleModelVa4_2"]
      p2 <- maxt/(1+exp(par["CasaleModelVa0_2"]*(par["CasaleModelVa1_2"]-t)))
    }
    if (any(names(par) == "CasaleModelVa0")) {
      maxt <- (0.01-par["CasaleModelVa4"])/(1+exp(par["CasaleModelVa2"]*(par["CasaleModelVa3"]-t))) + par["CasaleModelVa4"]
      p1 <- p2 <- maxt/(1+exp(par["CasaleModelVa0"]*(par["CasaleModelVa1"]-t)))
    }
    
    # modèle Rivalan 2005
    if (any(names(par) == "a0_1")) {
      delta <- par["delta_1"]
      if (is.na(delta)) delta <- 0
      mint <- par["a4_1"]/(1+exp(par["a2_1"]*(par["a3_1"]-t)))
      p1 <- (1-mint)/(1+exp(par["a0_1"]*(par["a1_1"]-t))) + mint + delta
    }
    if (any(names(par) == "a0_2")) {
      delta <- par["delta_2"]
      if (is.na(delta)) delta <- 0
      mint <- par["a4_2"]/(1+exp(par["a2_2"]*(par["a3_2"]-t)))
      p2 <- (1-mint)/(1+exp(par["a0_2"]*(par["a1_2"]-t))) + mint + delta
    }
    if (any(names(par) == "a0")) {
      delta <- par["delta"]
      if (is.na(delta)) delta <- 0
      mint <- par["a4"]/(1+exp(par["a2"]*(par["a3"]-t)))
      p1 <- p2 <- (1-mint)/(1+exp(par["a0"]*(par["a1"]-t))) + mint + delta
    }
    
    suffixe <- c("", "_1", "_2", "_L1", "_R1", "_L2", "_R2")
    
    for (suf in suffixe) {
      
      if (any(names(par) == paste0("D1", suf))) {
        if (!is.na(par[paste0("D1", suf)])) D1 <- abs(par[paste0("D1", suf)]) else D1 <- 0
        if (!is.na(par[paste0("D2D1", suf)])) D2D1 <- abs(par[paste0("D2D1", suf)]) else D2D1 <- max(t)+1
        D2 <- D1 + D2D1
        if (!is.na(par[paste0("D3D2", suf)])) D3D2 <- abs(par[paste0("D3D2", suf)]) else D3D2 <- max(t)+1
        D3 <- D2 + D3D2
        if (!is.na(par[paste0("A", suf)])) A <- invlogit(-par[paste0("A", suf)]) else A <- 0
        if (A == 0) A <- 1E-9
        if (A == 1) A <- 1-1E-9
        if (!is.na(par[paste0("B", suf)])) B <- invlogit(-par[paste0("B", suf)]) else B <- 0
        if (B == 0) B <- 1E-9
        if (B == 1) B <- 1-1E-9
        if (!is.na(par[paste0("C", suf)])) C <- invlogit(-par[paste0("C", suf)]) else C <- 0
        if (C == 0) C <- 1E-9
        if (C == 1) C <- 1-1E-9
        
        if (is.na(par[paste0("delta", suf)])) delta <- 0 else delta <- par[paste0("delta", suf)]
        if (D1==0) D1_2 <- 0.1 else D1_2 <- 2*D1
        if (D3==D2) D3D2 <- (D3-0.1)*D2 else D3D2 <-D3-D2
        
        suffi <- gsub("_", "", suf)
        if (suf=="") {
          p1 <- p2 <- ifelse(t<=D1, 
                             (1+cos(pi*((t+D1)/(D1_2))))*(A-B)+B, 
                             ifelse(t<=D2, B, 
                                    ifelse(t<=D3, (1+cos(pi*((t-D2)/(D3D2))))*0.5*(B-C)+C, C))) + delta
        } else {
          assign(paste0("p", suffi),  ifelse(t<=D1, 
                                             (1+cos(pi*((t+D1)/(D1_2))))*(A-B)+B, 
                                             ifelse(t<=D2, B, 
                                                    ifelse(t<=D3, (1+cos(pi*((t-D2)/(D3D2))))*0.5*(B-C)+C, C))) + delta)
        }
      }
    }
    
    if (!is.null(model_after)) eval(parse(text=model_after), envir= environment())
    
    if (is.null(model)) {
      z <- p1
    } else {
      z <- get(paste0("p", model))
    }
    
    return(as.numeric(ifelse(z<0, 1E-9, ifelse(z>1, 1-1E-9, z))))
    
    
  }
  
  # Si je suis ici c'est que je dois calculer une erreur standard
  
  if (method == "delta") {
    if (!(is.element('nlWaldTest', installed.packages()[,1]))) stop("nlWaldTest package must be installed to use the delata method.")
    
    suppressPackageStartupMessages(requireNamespace("nlWaldTest"))
    message("Estimation of distribution using delta method")
    
    VCov <- solve(Hessian)
    par_hess <- par
    par_hess <- par_hess[colnames(VCov)]
    
    par_add <- par
    par_add <- par_add[!names(par_add) %in% colnames(VCov)]
    if (identical(par_add,structure(numeric(0), .Names = character(0)))) par_add <- NULL
    
    gh <- data.frame(time=numeric(), value=numeric(), 
                     mean=numeric(), se=numeric(), 
                     'X2.5'=numeric(), 'X97.5'=numeric())
    
    if (any(diag(VCov)<0)) stop("Hessian matrix is not correct; probably you are not at maximum likelihood.")
    
    
    try_g2 <- function(..., Time, model_before=NULL, model_after=NULL, model, pfixed=NULL) {
      
      par <- c(...)
      
      # print(d(par))
      
      tgm <- Tagloss_model(Time, par=c(par, pfixed), 
                           model_before=model_before,
                           model_after=model_after, 
                           Hessian = NULL, model = model)
      
      return(tgm)
    }
    
    
    nlConfint2 <- function (obj = NULL, texts, level = 0.95, coeff = NULL, Vcov = NULL, 
                            df2 = NULL, x = NULL, parameters=NULL) 
    {
      if (!is.null(obj)) {
        co = try(coef(obj), silent = T)
        cond = attr(co, "condition")
        if (is.null(coeff) && (is.null(cond))) 
          coeff = co
        vc = try(vcov(obj), silent = T)
        cond2 = attr(vc, "condition")
        if (is.null(Vcov) && (is.null(cond2))) 
          Vcov = vc
      }
      if (is.null(coeff)) {
        if (is.null(obj)) 
          mess = "Both  'obj' and 'coeff' are missing"
        else {
          clm <- class(obj)
          part1 <- "There are no coef() methods for model objects of class \""
          mess <- paste0(part1, clm, "\".\nInput the 'coeff' parameter.")
        }
        stop(mess)
      }
      if (is.null(Vcov)) {
        if (is.null(obj)) 
          mess = "Both  'obj' and 'Vcov' are missing"
        else {
          clm <- class(obj)
          part1 <- "There are no vcov() methods for model objects of class \""
          mess <- paste0(part1, clm, "\".\nInput the 'Vcov' parameter.")
        }
        stop(mess)
      }
      if (length(texts) > 1) { kkk = texts[1]
      } else kkk = strsplit(texts[1], ";")[[1]]
      kkkfl = as.formula(paste("~", kkk[1]))
      vvss = setdiff(all.vars(kkkfl), "x")
      texts = getFromNamespace(".smartsub", ns="nlWaldTest")(vvss, "b", texts)
      if (length(texts) > 1) { ltext0 = texts
      } else ltext0 = strsplit(texts, ";")[[1]]
      texts1 = gsub("[", "", texts, fixed = T)
      texts1 = gsub("]", "", texts1, fixed = T)
      if (length(texts1) > 1) { ltext = texts1
      } else ltext = strsplit(texts1, ";")[[1]]
      r = length(ltext)
      n = length(coeff)
      namess = paste0("b", 1:n)
      for (j in 1L:n) assign(namess[j], coeff[j])
      if (!is.null(x)) {
        nx = length(x)
        namesx = paste0("x", 1:nx)
        for (j in 1L:nx) assign(namesx[j], x[j])
      }
      grad = c()
      hess = c()
      for (i in 1L:r) {
        if (!is.null(parameters)) {
          fli <- as.formula(paste("~", paste0(gsub(")", "", ltext[i]), ", ", 
                                              parameters, ")")))
        } else {
          fli <- as.formula(paste("~", ltext[i]))
        }
        z = try(deriv(as.formula(fli), namess), silent = T)
        if (inherits(z, "try-error")) {
          tei = as.character(i)
          tri2 = ", numerical derivatives were used in delta-method"
          wate = paste0("Note: For function ", i, tri2)
          message(wate)
          if (!is.null(parameters)) {
            ez = numericDeriv(quote(eval(parse(text = paste0(gsub(")", "", ltext[i]), ", ", 
                                                             parameters, ")")))), 
                              theta=namess)
          } else {
            ez = numericDeriv(quote(eval(parse(text = ltext[i]))), 
                              theta=namess)
          }
        } else ez = eval(z)
        hessj = attr(ez, "gradient")
        grad = rbind(grad, ez[1])
        hess = rbind(hess, hessj)
      }
      Rb = grad
      
      
      ddd = hess %*% Vcov %*% t(hess)
      # On retire ce test
      # matr = chol2inv(chol(ddd))
      ses = sqrt(diag(ddd))
      
      
      
      trydf = identical(df2, T)
      if (trydf) {
        isdf = try(df.residual(obj), silent = T)
        df2 = isdf
        if (is.null(df2)) {
          wn = "Note: Failed to extract df for denomenator; z-intervals applied"
          message(wn)
        }
      }
      getFromNamespace(".getint", ns="nlWaldTest")(as.numeric(Rb), ltext0, ses, level, df = df2)
    }
    
    for (ti in t) {  
      
      if (is.null(model_before)) {
        mb <- "model_before=NULL"
      } else {
        mb <- paste0("model_before=\"", model_before, "\"")
      }
      
      if (is.null(model_after)) {
        mb <- paste0(mb, ", model_after=NULL")
      } else {
        mb <- paste0(mb, ", model_after=\"", model_after, "\"")
      }
      
      if (!is.null(model)) {
        mb <- paste0(mb, ", model=\"", model, "\"")
      }
      
      if (!is.null(par_add)) {
        mb <- paste0(mb, ", pfixed=")
        mb <- paste0(mb, paste0("c(", paste(names(par_add), "=", unname(par_add), collapse = ", "), ")"))
      }
      
      
      ghi <- suppressMessages ( nlConfint2(texts=paste0(c("try_g2(", paste0("b[", seq_along(par_hess), "], ", collapse=""),"Time=", ti, ")"), 
                                                        collapse = ""), 
                                           parameters = mb, 
                                           level = 0.95, coeff = par_hess,
                                           Vcov = VCov, df2 = TRUE)
      )
      
      gh <- rbind(gh, data.frame(time=ti, value=ghi[, 1], mean=NA, sd=NA, 
                                 'X2.5'=ghi[, 2], 'X97.5'=ghi[, 3]))
      
    }
    
    return(gh)
  }
  
  # Ce n'est pas null et pas delta; je peux utiliser RandomFromHessianOrMCMC
  
  dfr <- RandomFromHessianOrMCMC(method = method, mcmc=mcmc, Hessian = Hessian, 
                                 fitted.parameters = par, replicates = replicates, 
                                 fn=Tagloss_model, ParTofn = "par",  model_before = model_before, 
                                 model_after = model_after, 
                                 model=model, x=x, t=t)
  return(dfr$quantiles)
}
