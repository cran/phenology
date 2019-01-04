#' Tagloss_cumul returns the cumulative rate of tag loss.
#' @title Return the cumulative rate of tag loss.
#' @author Marc Girondot
#' @return Return the cumulative rate of tag loss if hessian is null or a data.frame with distribution of cumulative rate of tag loss if hessian is not null.
#' @param t Time for which values of model must be estimated
#' @param x A Tagloss fitted model
#' @param par Parameters
#' @param hessian hessian matrix of parameters
#' @param model_before Function to be used before estimation of daily tagloss rate
#' @param model_after Function to be used after estimation of daily tagloss rate
#' @param model The model of parameter to use, can be N2, N1 or N0; or NLR, NL0, N0R, or N00 or NULL if hessian is NULL.
#' @param replicates Number of replicates to estimate se of output for resampling method
#' @description This function compute a model of cumulative tag loss rate for days t 
#' based on a set of parameters, par.\cr
#' If hessian is not null, it will estimate standard error of the output using numerical delta method is replicates 
#' is null or using resampling if replicates is not null.\cr
#' Parameters are described in \code{\link{Tagloss_fit}}.
#' @family Model of Tag-loss
#' @examples
#' \dontrun{
#' # Example
#' library(phenology)
#' 
#' # Data from Rivalan et al. 2005 - Table 2, line 1 - Fig 1D
#' par <- c(a0_2=-5.43E-2, a1_2=-103.52, a2_2=0, a3_2=0, a4_2=5.62E-4)
#' (y <- Tagloss_cumul(t=(1:6)*365, par=par))
#' plot(y[, "time"], y[, "N2"], type="l", bty="n", 
#'      xlab="Days after tagging", ylab="N2 proportion")
#' 
#' # Data from Rivalan et al. 2005 - Table 2, line 2 - Fig 1E
#' par <- c(a0_2=-6.80E-2, a1_2=-81.15, a2_2=-2.20E-4, a3_2=6348.01, a4_2=1.65E-3)
#' (y <- Tagloss_cumul(t=(1:6)*365, par=par))
#' plot(y[, "time"], y[, "N2"], type="l", 
#'      xlab="Days after tagging", ylab="N2 proportion")
#' 
#' # Data from Rivalan et al. 2005 - Table 2, line 3 - Fig 1F
#' par <- c(a0_2=-6.93E-2, a1_2=-78.92, a2_2=8.45E-4, a3_2=-16272.76, a4_2=2.87E-4)
#' (y <- Tagloss_cumul(t=(1:6)*365, par=par))
#' plot(y[, "time"], y[, "N2"], type="l", 
#'      xlab="Days after tagging", ylab="N2 proportion")
#' 
#' # Data from Rivalan et al. 2005 - Table 2, line 4 - Fig 1C
#' par <- c(a0_2=-1.68E-3, a1_2=-4141.68, a2_2=0, a3_2=0, a4_2=0)
#' (y <- Tagloss_cumul(t=(1:6)*365, par=par))
#' plot(y[, "time"], y[, "N2"], type="l", 
#'      xlab="Days after tagging", ylab="N2 proportion")
#' 
#' # Data from Rivalan et al. 2005 - Table 2, line 5 - Fig 1B
#' par <- c(a0_2=-3.77E-4, a1_2=-2000, a2_2=-0.001, a3_2=0, a4_2=5.00E-8)
#' (y <- Tagloss_cumul(t=(1:6)*365, par=par))
#' plot(y[, "time"], y[, "N2"], type="l", 
#'      xlab="Days after tagging", ylab="N2 proportion")
#' 
#' # Data from Rivalan et al. 2005 - Table 2, line 6 - Fig 1A
#' par <- c(a0_2=-1E5, a1_2=-2000, a2_2=0, a3_2=4000, a4_2=8.34E-4)
#' (y <- Tagloss_cumul(t=(1:6)*365, par=par))
#' plot(y[, "time"], y[, "N2"], type="l", 
#'      xlab="Days after tagging", ylab="N2 proportion")
#'      
#' # Data from Rivalan et al. 2005 - Table 2, line 1 - Fig 1D
#' # With tagloss rate dependency on tage number
#' par <- c(a0_2=-5.43E-2, a1_2=-103.52, a2_2=0, a3_2=0, a4_2=5.62E-4, 
#'          a0_1=-5.43E-2, a1_1=-103.52, a2_1=0, a3_1=0, a4_1=5.62E-4, delta_1=3.2E-4)
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, 
#'                          model="Cumul", 
#'                          decoration = TRUE)
#' 
#' p2 <- Tagloss_model(t=1:(6*365), par=par, model="2")
#' p1 <- Tagloss_model(t=1:(6*365), par=par, model="1")
#' par(mar=c(4, 5, 2, 1))
#' plot(x=1:(6*365), y=p2, bty="n", type="l", las=1, ylim=c(0,0.003), ylab="")
#' mtext("Daily tag loss", side=2, line=4)
#' lines(x=1:(6*365), y=p1, col="red")
#' legend("topright", legend=c("2>1", "1>0"), lty=1, col=c("black", "red"))
#' 
#' Tagloss_cumul(t=(1:6)*365, par=par)
#' 
#' # Without tagloss rate dependency on tag number
#' par <- c(a0_2=-5.43E-2, a1_2=-103.52, a2_2=0, a3_2=0, a4_2=5.62E-4, 
#'          a0_1=-5.43E-2, a1_1=-103.52, a2_1=0, a3_1=0, a4_1=5.62E-4)
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, 
#'                          model="Cumul", 
#'                          decoration = TRUE)
#' Tagloss_cumul(t=(1:6)*365, par=par)
#' 
#' #### Data from Casale et al. 2017
#' # Table 1 - Model II
#' par <- c(CasaleModelIIa0_2=-0.0511, CasaleModelIIa1_2=-100, CasaleModelIIa4_2=0.00014)
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, 
#'                          model="Cumul", 
#'                          decoration = TRUE)
#' Tagloss_cumul(t=(1:6)*365, par=par)
#' 
#' # Table 1 - Model IV
#' par <- c(CasaleModelIVa0_2=-0.0132, CasaleModelIVa1_2=-100, 
#'          CasaleModelIVa2_2=0.0327, CasaleModelIVa3_2=109.98, 
#'          CasaleModelIVa4_2=0.00011)
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, 
#'                          model="Cumul", 
#'                          decoration = TRUE)
#' Tagloss_cumul(t=(1:6)*365, par=par)
#' 
#' # Table 1 - Model I
#' par <- c(CasaleModelIc_2=0.00027)
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, 
#'                          model="Cumul", 
#'                          decoration = TRUE)
#' Tagloss_cumul(t=(1:6)*365, par=par)
#' 
#' # Table 1 - Model III
#' par <- c(CasaleModelIIIa0_2=1.14E-10, CasaleModelIIIa1_2=-110.04, 
#'          CasaleModelIIIa4_2=0.00055)
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, 
#'                          model="Cumul", 
#'                          decoration = TRUE)
#' Tagloss_cumul(t=(1:6)*365, par=par)
#' 
#' # Table 1 - Model V
#' par <- c(CasaleModelVa0_2=4.04E-10, CasaleModelVa1_2=-90, 
#'          CasaleModelVa2_2=-0.0326, CasaleModelVa3_2=100.31, 
#'          CasaleModelVa4_2=0.00006)
#' phenology:::plot.Tagloss(x=list(), t=1:1000, fitted.parameters=par, 
#'                          model="Cumul", 
#'                          decoration = TRUE)
#' Tagloss_cumul(t=(1:6)*365, par=par)
#' 
#' }
#' @export

Tagloss_cumul <- function(t, par=NULL, hessian=NULL, model_before = NULL, 
                          model_after=NULL, 
                          model=NULL, replicates=NULL, x=NULL) {
  
  # par=NULL; hessian=NULL; model_before = NULL; model_after=NULL; model="cumul"; replicates=NULL; x=NULL
  
  if (!is.null(model)) model <- toupper(model)
  # D'abord je fais sans hessian
  
  
  
  if (is.null(x) & is.null(par)) {
    stop("Both par and x cannot be null at the same time")
  }
  
  if (!is.null(x)) {
    if (is.null(par)) par <- c(x$par, x$fixed.par)
    if (is.null(hessian)) hessian <- x$hessian
    if (is.null(model_before)) model_before <- x$model_before
    if (is.null(model_after)) model_after <- x$model_after
  }
  
  days.maximum <- max(t)
  
  if (is.null(hessian)) {
    
    if (!is.null(model_before)) eval(parse(text=model_before), envir= environment())
    
    if (any(grepl("_R2", names(par)))) pR2 <- Tagloss_model(1:days.maximum, par, model_before = model_before, model_after = model_after, model="R2") else pR2 <- NA
    if (any(grepl("_R1", names(par)))) pR1 <- Tagloss_model(1:days.maximum, par, model_before = model_before, model_after = model_after, model="R1") else pR1 <- NA
    if (any(grepl("_L2", names(par)))) pL2 <- Tagloss_model(1:days.maximum, par, model_before = model_before, model_after = model_after, model="L2") else pL2 <- NA
    if (any(grepl("_L1", names(par)))) pL1 <- Tagloss_model(1:days.maximum, par, model_before = model_before, model_after = model_after, model="L1") else pL1 <- NA
    
    if (any(grepl("_2", names(par)))) p2 <- Tagloss_model(1:days.maximum, par, model_before = model_before, model_after = model_after, model="2") else p2 <- NA
    if (any(grepl("_1", names(par)))) p1 <- Tagloss_model(1:days.maximum, par, model_before = model_before, model_after = model_after, model="1") else p1 <- NA
    
    if (!any(grepl("_", names(par)))) p1 <- p2 <- Tagloss_model(1:days.maximum, par, model_before = model_before, model_after = model_after, model=NULL)
    
    if (!is.null(model_after)) eval(parse(text=model_after), envir= environment())
    
    if (is.na(p1[1]) & !is.na(p2[1])) p1 <- p2
    if (is.na(p2[1]) & !is.na(p1[1])) p2 <- p1
    
    if (is.na(pR1[1]) & !is.na(pR2[1])) pR1 <- pR2
    if (is.na(pR2[1]) & !is.na(pR1[1])) pR2 <- pR1
    
    if (is.na(pL1[1]) & !is.na(pL2[1])) pL1 <- pL2
    if (is.na(pL2[1]) & !is.na(pL1[1])) pL2 <- pL1
    
    
    if ((!is.na(p1[1])) & (!is.na(p2[1]))) {
      ################
      # modèle sans RL
      ################
      
      d <- data.frame(time=1:days.maximum, N2=rep(x = 1, times=days.maximum), N1=rep(x=0, times=days.maximum), N0=rep(x=0, times=days.maximum))
      if (days.maximum >1) {
        for (i in 2:days.maximum) {
          d[i, "N2"] <- d[i-1, "N2"] * (1-p2[i-1])
          d[i, "N1"] <- d[i-1, "N1"] * (1-p1[i-1]) + d[i-1, "N2"] * p2[i-1] * (1-p1[i-1])
          d[i, "N0"] <- d[i-1, "N0"] + d[i-1, "N1"] * p1[i-1] + d[i-1, "N2"] * p2[i-1] * p1[i-1]
        }
      }
      d[, -1] <- d[, -1]/rowSums(d[, -1])
      
      
    } else {
      ################
      # modèle avec RL
      ################
      
      d <- data.frame(time=1:days.maximum, NLR=rep(x = 1, times=days.maximum), NL0=rep(x=0, times=days.maximum), N0R=rep(x=0, times=days.maximum), N00=rep(x=0, times=days.maximum))
      if (days.maximum >1) {
        for (i in 2:days.maximum) {
          d[i, "NLR"] <- d[i-1, "NLR"] * (1-pR2[i-1]) * (1-pL2[i-1])
          d[i, "NL0"] <- d[i-1, "NL0"] * (1-pR1[i-1]) + d[i-1, "NLR"] * pL2[i-1] * (1-pR1[i-1])
          d[i, "N0R"] <- d[i-1, "N0R"] * (1-pL1[i-1]) + d[i-1, "NLR"] * pR2[i-1] * (1-pL1[i-1])
          d[i, "N00"] <- d[i-1, "N00"] + d[i-1, "NL0"] * pL1[i-1] + d[i-1, "N0R"] * pR1[i-1] + d[i-1, "NLR"] * pL2[i-1] * pR1[i-1] + d[i-1, "NLR"] * pR2[i-1] * pL1[i-1]
        }
      }
      d[, -1] <- d[, -1]/rowSums(d[, -1])
    }
    
    d <- d[match(t, d[, "time"]), ]
    return(d)
    
    
  } else {
    
    # Bon, j'ai une matrice... mais maintenant je dois l'utliser !
    
    VCov <- solve(hessian)
    
    par_hess <- par[colnames(VCov)]
    par_add <- par[!names(par) %in% colnames(VCov)]
    if (identical(par_add,structure(numeric(0), .Names = character(0)))) par_add <- NULL
    
    gh <- data.frame(time=numeric(), value=numeric(), 
                     mean=numeric(), se=numeric(), 
                     'X2.5'=numeric(), 'X97.5'=numeric())
    
    
    suppressPackageStartupMessages(  requireNamespace("nlWaldTest"))
    message("Estimation of distribution using delta method")
    
    try_N <- function(..., Time, model_before=NULL, model_after=NULL, pfixed=NULL, model) {
      
      par <- c(...)
      
      tgm <- Tagloss_cumul(Time, par=c(par, pfixed), 
                           model_before=model_before, 
                           model_after=model_after, 
                           hessian = NULL, model = NULL)
      
      return(tgm[1, model])
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
          clm = class(obj)
          part1 = "There are no coef() methods for model objects of class \""
          mess = paste0(part1, clm, "\".\nInput the 'coeff' parameter.")
        }
        stop(mess)
      }
      if (is.null(Vcov)) {
        if (is.null(obj)) 
          mess = "Both  'obj' and 'Vcov' are missing"
        else {
          clm = class(obj)
          part1 = "There are no vcov() methods for model objects of class \""
          mess = paste0(part1, clm, "\".\nInput the 'Vcov' parameter.")
        }
        stop(mess)
      }
      if (length(texts) > 1) 
      { kkk = texts[1]
      } else kkk = strsplit(texts[1], ";")[[1]]
      kkkfl = as.formula(paste("~", kkk[1]))
      vvss = setdiff(all.vars(kkkfl), "x")
      texts = getFromNamespace(".smartsub", ns="nlWaldTest")(vvss, "b", texts)
      if (length(texts) > 1) 
      { ltext0 = texts
      } else ltext0 = strsplit(texts, ";")[[1]]
      texts1 = gsub("[", "", texts, fixed = T)
      texts1 = gsub("]", "", texts1, fixed = T)
      if (length(texts1) > 1) 
      { ltext = texts1
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
        if (class(z) == "try-error") {
          tei = as.character(i)
          tri2 = ", numerical derivatives were used in delta-method"
          wate = paste0("Note: For function ", i, tri2)
          message(wate)
          if (!is.null(parameters)) {
            ez = numericDeriv(quote(eval(parse(text = paste0(gsub(")", "", ltext[i]), ", ", 
                                                             parameters, ")")))), 
                              namess)
          } else {
            ez = numericDeriv(quote(eval(parse(text = ltext[i]))), 
                              namess)
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
      
      mb <- paste0(mb, ", model=\"", model, "\"")
      
      if (!is.null(par_add)) {
        mb <- paste0(mb, ", pfixed=")
        mb <- paste0(mb, paste0("c(", paste(names(par_add), "=", unname(par_add), collapse = ", "), ")"))
      }
      
      
      ghi <- suppressMessages ( nlConfint2(texts=paste0(c("try_N(", paste0("b[", seq_along(par_hess), "], ", collapse=""),"Time=", ti, ")"), 
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
  
}

