#' AutoFitPhenology runs fit for phenology and tests several combinations
#' @title Automatic fit for phenology and tests
#' @author Marc Girondot
#' @return A list with 12 elements corresponding to the 12 tested models
#' @param data Dataset generated with add_phenology()
#' @param progressbar If FALSE, do not show the progress bar
#' @param ... Parameters for fit_phenology()
#' @description This function is used to test several combinations of fit at a time.
#' @family Phenology model
#' @examples
#' \dontrun{
#' library(phenology)
#' # Read a file with data
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot <- add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Run the optimisation
#' result_Gratiot_Auto <- AutoFitPhenology(data=data_Gratiot)
#' result_Gratiot_Auto <- AutoFitPhenology(data=data_Gratiot, 
#'           control=list(trace=0, REPORT=100, maxit=500))
#' }
#' @export


AutoFitPhenology <-
  function(data=stop("A dataset must be provided"), 
           progressbar=TRUE, 
           ...) {
    
    if (!inherits(data, "phenologydata")) {
      stop("Data must be formated first using the function add_phenology().")
    }
    
    p3p <- list(...)
    
    if (progressbar) pb<-txtProgressBar(min=1, max=8, style=3)
    
################ LengthB, LengthE
    
    if (progressbar) setTxtProgressBar(pb, 1)
    # LB, LE, PK, FL, PMinB, PMinE
    pfixed <- c(Min=0)
    parg <- par_init(data = data, fixed.parameters = pfixed)
    parg <- c(parg, PMinB=5, PMinE=5)
    pfixed <- NULL
    
    p3p_encours <- modifyList(p3p, list(data = data, 
                                        fitted.parameters = parg, 
                                        fixed.parameters = pfixed,
                                        hessian = FALSE))
    fit1 <- do.call("fit_phenology", p3p_encours)
    
    if (progressbar) setTxtProgressBar(pb, 2)
    # LB, LE, PK, FL, PMin
    pfixed <- NULL
    parg <- fit1$par
    parg["PMinB"] <- mean(c(parg["PMinB"], parg["PMinE"]))
    parg <- parg[-which(names(parg)=="PMinE")]
    names(parg)[which(names(parg)=="PMinB")] <- "PMin"
    
    p3p_encours <- modifyList(p3p, list(data = data, 
                                        fitted.parameters = parg, 
                                        fixed.parameters = pfixed,
                                        hessian = FALSE))
    fit2 <- do.call("fit_phenology", p3p_encours)
    
    if (progressbar) setTxtProgressBar(pb, 3)
    # LB, LE, PK, FL
    pfixed <- c(Min = 0)
    parg <- fit2$par
    parg <- parg[-which(names(parg)=="PMin")]

    p3p_encours <- modifyList(p3p, list(data = data, 
                                        fitted.parameters = parg, 
                                        fixed.parameters = pfixed,
                                        hessian = FALSE))
    fit3 <- do.call("fit_phenology", p3p_encours)
    
    if (progressbar) setTxtProgressBar(pb, 4)
    # LB, LE, PK, PMinB, PMinE
    pfixed <- c(Flat = 0)
    parg <- fit1$par
    parg <- parg[-which(names(parg)=="Flat")]
    
    p3p_encours <- modifyList(p3p, list(data = data, 
                                        fitted.parameters = parg, 
                                        fixed.parameters = pfixed,
                                        hessian = FALSE))
    fit4 <- do.call("fit_phenology", p3p_encours)
    
    if (progressbar) setTxtProgressBar(pb, 5)
    # LB, LE, PK, PMin
    pfixed <- c(Flat = 0)
    parg <- fit4$par
    parg["PMinB"] <- mean(c(parg["PMinB"], parg["PMinE"]))
    parg <- parg[-which(names(parg)=="PMinE")]
    names(parg)[which(names(parg)=="PMinB")] <- "PMin"
    
    p3p_encours <- modifyList(p3p, list(data = data, 
                                        fitted.parameters = parg, 
                                        fixed.parameters = pfixed,
                                        hessian = FALSE))
    fit5 <- do.call("fit_phenology", p3p_encours)

    if (progressbar) setTxtProgressBar(pb, 6)
    # LB, LE, PK
    pfixed <- c(Min = 0, Flat = 0)
    parg <- fit5$par
    parg <- parg[-which(names(parg)=="PMin")]
    
    p3p_encours <- modifyList(p3p, list(data = data, 
                                        fitted.parameters = parg, 
                                        fixed.parameters = pfixed,
                                        hessian = FALSE))
    fit6 <- do.call("fit_phenology", p3p_encours)
    
############## L seulement
    if (progressbar) setTxtProgressBar(pb, 7)
    # LB, PK, FL, PMinB, PMinE
    parg <- fit1$par
    parg["LengthB"] <- mean(c(parg["LengthB"], parg["LengthE"]))
    parg <- parg[-which(names(parg)=="LengthE")]
    names(parg)[which(names(parg)=="LengthB")] <- "Length"
    pfixed <- fit1$fixed.parameters
    
    p3p_encours <- modifyList(p3p, list(data = data, 
                                        fitted.parameters = parg, 
                                        fixed.parameters = pfixed,
                                        hessian = FALSE))
    fit7 <- do.call("fit_phenology", p3p_encours)

    if (progressbar) setTxtProgressBar(pb, 8)
    # L, PK, FL, PMin
    parg <- fit2$par
    parg["LengthB"] <- mean(c(parg["LengthB"], parg["LengthE"]))
    parg <- parg[-which(names(parg)=="LengthE")]
    names(parg)[which(names(parg)=="LengthB")] <- "Length"
    pfixed <- fit2$fixed.parameters
    
    p3p_encours <- modifyList(p3p, list(data = data, 
                                        fitted.parameters = parg, 
                                        fixed.parameters = pfixed,
                                        hessian = FALSE))
    fit8 <- do.call("fit_phenology", p3p_encours)
    
    if (progressbar) setTxtProgressBar(pb, 9)
    # L, PK, FL
    parg <- fit3$par
    parg["LengthB"] <- mean(c(parg["LengthB"], parg["LengthE"]))
    parg <- parg[-which(names(parg)=="LengthE")]
    names(parg)[which(names(parg)=="LengthB")] <- "Length"
    pfixed <- fit3$fixed.parameters
    
    p3p_encours <- modifyList(p3p, list(data = data, 
                                        fitted.parameters = parg, 
                                        fixed.parameters = pfixed,
                                        hessian = FALSE))
    fit9 <- do.call("fit_phenology", p3p_encours)
    
    if (progressbar) setTxtProgressBar(pb, 10)
    # L, PK, PMinB, PMinE
    parg <- fit4$par
    parg["LengthB"] <- mean(c(parg["LengthB"], parg["LengthE"]))
    parg <- parg[-which(names(parg)=="LengthE")]
    names(parg)[which(names(parg)=="LengthB")] <- "Length"
    pfixed <- fit4$fixed.parameters
    
    p3p_encours <- modifyList(p3p, list(data = data, 
                                        fitted.parameters = parg, 
                                        fixed.parameters = pfixed,
                                        hessian = FALSE))
    fit10 <- do.call("fit_phenology", p3p_encours)
    
    if (progressbar) setTxtProgressBar(pb, 11)
    # L, PK, PMin
    parg <- fit5$par
    parg["LengthB"] <- mean(c(parg["LengthB"], parg["LengthE"]))
    parg <- parg[-which(names(parg)=="LengthE")]
    names(parg)[which(names(parg)=="LengthB")] <- "Length"
    pfixed <- fit5$fixed.parameters
    
    p3p_encours <- modifyList(p3p, list(data = data, 
                                        fitted.parameters = parg, 
                                        fixed.parameters = pfixed,
                                        hessian = FALSE))
    fit11 <- do.call("fit_phenology", p3p_encours)
    
    if (progressbar) setTxtProgressBar(pb, 12)
    # LB, LE, PK, FL
    parg <- fit6$par
    parg["LengthB"] <- mean(c(parg["LengthB"], parg["LengthE"]))
    parg <- parg[-which(names(parg)=="LengthE")]
    names(parg)[which(names(parg)=="LengthB")] <- "Length"
    pfixed <- fit6$fixed.parameters
    
    p3p_encours <- modifyList(p3p, list(data = data, 
                                        fitted.parameters = parg, 
                                        fixed.parameters = pfixed,
                                        hessian = FALSE))
    fit12 <- do.call("fit_phenology", p3p_encours)
    
    compare_AIC(LengthB_LengthE_Flat_PMinB_PMinE = fit1, 
                LengthB_LengthE_Flat_PMin = fit2, 
                LengthB_LengthE_Flat = fit3, 
                LengthB_LengthE_PMinB_PMinE = fit4, 
                LengthB_LengthE_PMin = fit5, 
                LengthB_LengthE = fit6, 
                Length_Flat_PMinB_PMinE = fit7, 
                Length_Flat_PMin = fit8, 
                Length_Flat = fit9, 
                Length_PMinB_PMinE = fit10, 
                Length_PMin = fit11, 
                Length = fit12)
    
    return(list(LengthB_LengthE_Flat_PMinB_PMinE = fit1, 
                LengthB_LengthE_Flat_PMin = fit2, 
                LengthB_LengthE_Flat = fit3, 
                LengthB_LengthE_PMinB_PMinE = fit4, 
                LengthB_LengthE_PMin = fit5, 
                LengthB_LengthE = fit6, 
                Length_Flat_PMinB_PMinE = fit7, 
                Length_Flat_PMin = fit8, 
                Length_Flat = fit9, 
                Length_PMinB_PMinE = fit10, 
                Length_PMin = fit11, 
                Length = fit12))
  }
