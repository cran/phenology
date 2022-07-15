#' Tagloss_mcmc_p generates set of parameters to be used with Tagloss_mcmc()
#' @title Generates set of parameters to be used with Tagloss_mcmc()
#' @author Marc Girondot
#' @return A matrix with the parameters
#' @param result An object obtained after a Tagloss_fit() fit
#' @param default.density The default density, "dnorm" or "dunif'
#' @param accept If TRUE, does not wait for use interaction
#' @description Interactive (or not!) script used to generate set of parameters to be used with Tagloss_mcmc().
#' @family Model of Tag-loss
#' @examples 
#' \dontrun{
#' library(phenology)
#' # Example
#' data_f_21 <- Tagloss_format(outLR, model="21")
#' 
#' # model fitted by Rivalan et al. 2005
#' par <- c(a0_2=-5.43E-2, a1_2=-103.52, a4_2=5.62E-4, 
#'          delta_1=3.2E-4)
#' pfixed <- c(a2_2=0, a3_2=0, a2_1=0, a3_1=0)
#' model_before <- "par['a0_1']=par['a0_2'];par['a1_1']=par['a1_2'];par['a4_1']=par['a4_2']"
#' o <- Tagloss_fit(data=data_f_21, fitted.parameters=par, fixed.parameters=pfixed, 
#'                  model_before=model_before)
#' pMCMC <- Tagloss_mcmc_p(o, accept=TRUE)
#' o_MCMC <- Tagloss_mcmc(data=data_f_21, parameters=pMCMC, fixed.parameters=pfixed, 
#'                  model_before=model_before, 
#'                  n.iter=10000, n.chains = 1, n.adapt = 100, thin=30)
#' }
#' @export

Tagloss_mcmc_p <- function(result=stop("An output from Tagloss_fit() must be provided"), 
                             default.density="dunif", 
                             accept=FALSE) {
  
  if (!inherits(result, "Tagloss")) {
    stop("An output from Tagloss_fit() must be provided")
  }
  
  # d'abord je sors les paramètres à utiliser
  
  par <- result$par
  
  # "D1_L2"
  pe <- ifelse(is.na(par["D1_L2"]), 50, par["D1_L2"])
  if (default.density == "dunif") {
    D1_L2 <-  c("dunif", 0, max(c(100, pe+2*pe)), 5, 0, max(c(100, pe+2*pe)), pe)
  } else {
    D1_L2 <-  c("dnorm", pe, pe/2, 5, 0, max(c(100, pe+2*pe)), pe)
  }
  # "D1_R2"
  pe <- ifelse(is.na(par["D1_R2"]), 50, par["D1_R2"])
  if (default.density == "dunif") {
    D1_R2 <-  c("dunif", 0, max(c(100, pe+2*pe)), 5, 0, max(c(100, pe+2*pe)), pe)
  } else {
    D1_R2 <-  c("dnorm", pe, pe/2, 5, 0, max(c(100, pe+2*pe)), pe)
  }
  # "D1_L1"
  pe <- ifelse(is.na(par["D1_L1"]), 50, par["D1_L1"])
  if (default.density == "dunif") {
    D1_L1 <-  c("dunif", 0, max(c(100, pe+2*pe)), 5, 0, max(c(100, pe+2*pe)), pe)
  } else {
    D1_L1 <-  c("dnorm", pe, pe/2, 5, 0, max(c(100, pe+2*pe)), pe)
  }
  # "D1_R1"
  pe <- ifelse(is.na(par["D1_R1"]), 50, par["D1_R1"])
  if (default.density == "dunif") {
    D1_R1 <-  c("dunif", 0, max(c(100, pe+2*pe)), 5, 0, max(c(100, pe+2*pe)), pe)
  } else {
    D1_R1 <-  c("dnorm", pe, pe/2, 5, 0, max(c(100, pe+2*pe)), pe)
  }
  # "D1_2"
  pe <- ifelse(is.na(par["D1_2"]), 50, par["D1_2"])
  if (default.density == "dunif") {
    D1_2 <-  c("dunif", 0, max(c(100, pe+2*pe)), 5, 0, max(c(100, pe+2*pe)), pe)
  } else {
    D1_2 <-  c("dnorm", pe, pe/2, 5, 0, max(c(100, pe+2*pe)), pe)
  }
  # "D1_1"
  pe <- ifelse(is.na(par["D1_1"]), 50, par["D1_1"])
  if (default.density == "dunif") {
    D1_1 <-  c("dunif", 0, max(c(100, pe+2*pe)), 5, 0, max(c(100, pe+2*pe)), pe)
  } else {
    D1_1 <-  c("dnorm", pe, pe/2, 5, 0, max(c(100, pe+2*pe)), pe)
  }
  
  # "D2D1_L2"
  pe <- ifelse(is.na(par["D2D1_L2"]), 200, par["D2D1_L2"])
  if (default.density == "dunif") {
    D2D1_L2 <-  c("dunif", 0, max(c(200, pe+2*pe)), 5, 0, max(c(200, pe+2*pe)), pe)
  } else {
    D2D1_L2 <-  c("dnorm", pe, pe/2, 10, 0, max(c(200, pe+2*pe)), pe)
  }
  # "D2D1_R2"
  pe <- ifelse(is.na(par["D2D1_R2"]), 200, par["D2D1_R2"])
  if (default.density == "dunif") {
    D2D1_R2 <-  c("dunif", 0, max(c(200, pe+2*pe)), 5, 0, max(c(200, pe+2*pe)), pe)
  } else {
    D2D1_R2 <-  c("dnorm", pe, pe/2, 10, 0, max(c(200, pe+2*pe)), pe)
  }
  # "D2D1_L1"
  pe <- ifelse(is.na(par["D2D1_L1"]), 200, par["D2D1_L1"])
  if (default.density == "dunif") {
    D2D1_L1 <-  c("dunif", 0, max(c(200, pe+2*pe)), 5, 0, max(c(200, pe+2*pe)), pe)
  } else {
    D2D1_L1 <-  c("dnorm", pe, pe/2, 10, 0, max(c(200, pe+2*pe)), pe)
  }
  # "D2D1_R1"
  pe <- ifelse(is.na(par["D2D1_R1"]), 200, par["D2D1_R1"])
  if (default.density == "dunif") {
    D2D1_R1 <-  c("dunif", 0, max(c(200, pe+2*pe)), 5, 0, max(c(200, pe+2*pe)), pe)
  } else {
    D2D1_R1 <-  c("dnorm", pe, pe/2, 10, 0, max(c(200, pe+2*pe)), pe)
  }
  # "D2D1_2"
  pe <- ifelse(is.na(par["D2D1_2"]), 200, par["D2D1_2"])
  if (default.density == "dunif") {
    D2D1_2 <-  c("dunif", 0, max(c(200, pe+2*pe)), 5, 0, max(c(200, pe+2*pe)), pe)
  } else {
    D2D1_2 <-  c("dnorm", pe, pe/2, 10, 0, max(c(200, pe+2*pe)), pe)
  }
  # "D2D1_1"
  pe <- ifelse(is.na(par["D2D1_1"]), 200, par["D2D1_1"])
  if (default.density == "dunif") {
    D2D1_1 <-  c("dunif", 0, max(c(200, pe+2*pe)), 5, 0, max(c(200, pe+2*pe)), pe)
  } else {
    D2D1_1 <-  c("dnorm", pe, pe/2, 10, 0, max(c(200, pe+2*pe)), pe)
  }
  
  # "D3D2_L2"
  pe <- ifelse(is.na(par["D3D2_L2"]), 200, par["D3D2_L2"])
  if (default.density == "dunif") {
    D3D2_L2 <-  c("dunif", 0, max(c(200, pe+2*pe)), 5, 0, max(c(200, pe+2*pe)), pe)
  } else {
    D3D2_L2 <-  c("dnorm", pe, pe/2, 10, 0, max(c(200, pe+2*pe)), pe)
  }
  # "D3D2_R2"
  pe <- ifelse(is.na(par["D3D2_R2"]), 200, par["D3D2_R2"])
  if (default.density == "dunif") {
    D3D2_R2 <-  c("dunif", 0, max(c(200, pe+2*pe)), 5, 0, max(c(200, pe+2*pe)), pe)
  } else {
    D3D2_R2 <-  c("dnorm", pe, pe/2, 10, 0, max(c(200, pe+2*pe)), pe)
  }
  # "D3D2_L1"
  pe <- ifelse(is.na(par["D3D2_L1"]), 200, par["D3D2_L1"])
  if (default.density == "dunif") {
    D3D2_L1 <-  c("dunif", 0, max(c(200, pe+2*pe)), 5, 0, max(c(200, pe+2*pe)), pe)
  } else {
    D3D2_L1 <-  c("dnorm", pe, pe/2, 10, 0, max(c(200, pe+2*pe)), pe)
  }
  # "D3D2_R1"
  pe <- ifelse(is.na(par["D3D2_R1"]), 200, par["D3D2_R1"])
  if (default.density == "dunif") {
    D3D2_R1 <-  c("dunif", 0, max(c(200, pe+2*pe)), 5, 0, max(c(200, pe+2*pe)), pe)
  } else {
    D3D2_R1 <-  c("dnorm", pe, pe/2, 10, 0, max(c(200, pe+2*pe)), pe)
  }
  # "D3D2_2"
  pe <- ifelse(is.na(par["D3D2_2"]), 200, par["D3D2_2"])
  if (default.density == "dunif") {
    D3D2_2 <-  c("dunif", 0, max(c(200, pe+2*pe)), 5, 0, max(c(200, pe+2*pe)), pe)
  } else {
    D3D2_2 <-  c("dnorm", pe, pe/2, 10, 0, max(c(200, pe+2*pe)), pe)
  }
  # "D3D2_1"
  pe <- ifelse(is.na(par["D3D2_1"]), 200, par["D3D2_1"])
  if (default.density == "dunif") {
    D3D2_1 <-  c("dunif", 0, max(c(200, pe+2*pe)), 5, 0, max(c(200, pe+2*pe)), pe)
  } else {
    D3D2_1 <-  c("dnorm", pe, pe/2, 10, 0, max(c(200, pe+2*pe)), pe)
  }
  
  
  # "A_L2"
  pe <- ifelse(is.na(par["A_L2"]), 6, par["A_L2"])
  if (default.density == "dunif") {
    A_L2 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    A_L2 <-  c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "A_R2"
  pe <- ifelse(is.na(par["A_R2"]), 6, par["A_R2"])
  if (default.density == "dunif") {
    A_R2 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    A_R2 <-  c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "A_L1"
  pe <- ifelse(is.na(par["A_L1"]), 6, par["A_L1"])
  if (default.density == "dunif") {
    A_L1 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    A_L1 <-  c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "A_R1"
  pe <- ifelse(is.na(par["A_R1"]), 6, par["A_R1"])
  if (default.density == "dunif") {
    A_R1 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    A_R1 <-  c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "A_2"
  pe <- ifelse(is.na(par["A_2"]), 6, par["A_2"])
  if (default.density == "dunif") {
    A_2 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    A_2 <-  c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "A_1"
  pe <- ifelse(is.na(par["A_1"]), 6, par["A_1"])
  if (default.density == "dunif") {
    A_1 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    A_1 <- c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  
  # "B_L2"
  pe <- ifelse(is.na(par["B_L2"]), 6, par["B_L2"])
  if (default.density == "dunif") {
    B_L2 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    B_L2 <-  c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "B_R2"
  pe <- ifelse(is.na(par["B_R2"]), 6, par["B_R2"])
  if (default.density == "dunif") {
    B_R2 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    B_R2 <-  c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "B_L1"
  pe <- ifelse(is.na(par["B_L1"]), 6, par["B_L1"])
  if (default.density == "dunif") {
    B_L1 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    B_L1 <-  c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "B_R1"
  pe <- ifelse(is.na(par["B_R1"]), 6, par["B_R1"])
  if (default.density == "dunif") {
    B_R1 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    B_R1 <-  c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "B_2"
  pe <- ifelse(is.na(par["B_2"]), 6, par["B_2"])
  if (default.density == "dunif") {
    B_2 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    B_2 <-  c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "B_1"
  pe <- ifelse(is.na(par["B_1"]), 6, par["B_1"])
  if (default.density == "dunif") {
    B_1 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    B_1 <- c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  
  # "C_L2"
  pe <- ifelse(is.na(par["C_L2"]), 6, par["C_L2"])
  if (default.density == "dunif") {
    C_L2 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    C_L2 <-  c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "C_R2"
  pe <- ifelse(is.na(par["C_R2"]), 6, par["C_R2"])
  if (default.density == "dunif") {
    C_R2 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    C_R2 <-  c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "C_L1"
  pe <- ifelse(is.na(par["C_L1"]), 6, par["C_L1"])
  if (default.density == "dunif") {
    C_L1 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    C_L1 <-  c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "C_R1"
  pe <- ifelse(is.na(par["C_R1"]), 6, par["C_R1"])
  if (default.density == "dunif") {
    C_R1 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    C_R1 <-  c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "C_2"
  pe <- ifelse(is.na(par["C_2"]), 6, par["C_2"])
  if (default.density == "dunif") {
    C_2 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    C_2 <-  c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "C_1"
  pe <- ifelse(is.na(par["C_1"]), 6, par["C_1"])
  if (default.density == "dunif") {
    C_1 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    C_1 <- c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  
  # "delta_L2"
  pe <- ifelse(is.na(par["delta_L2"]), 1, par["delta_L2"])
  if (default.density == "dunif") {
    delta_L2 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    delta_L2 <- c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "delta_R2"
  pe <- ifelse(is.na(par["delta_R2"]), 1, par["delta_R2"])
  if (default.density == "dunif") {
    delta_R2 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    delta_R2 <- c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "delta_L1"
  pe <- ifelse(is.na(par["delta_L1"]), 1, par["delta_L1"])
  if (default.density == "dunif") {
    delta_L1 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    delta_L1 <- c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "delta_R1"
  pe <- ifelse(is.na(par["delta_R1"]), 1, par["delta_R1"])
  if (default.density == "dunif") {
    delta_R1 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    delta_R1 <- c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "delta_2"
  pe <- ifelse(is.na(par["delta_2"]), 1, par["delta_2"])
  if (default.density == "dunif") {
    delta_2 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    delta_2 <- c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  # "delta_1"
  pe <- ifelse(is.na(par["delta_1"]), 1, par["delta_1"])
  if (default.density == "dunif") {
    delta_1 <-  c("dunif", min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  } else {
    delta_1 <- c("dnorm", pe, pe/2, 1, min(c(-20, pe-pe/10)), max(c(20, pe+pe/10)), pe)
  }
  
  priors <- list(D1_L2, D2D1_L2, D3D2_L2, A_L2, B_L2, C_L2, delta_L2,
                 D1_R2, D2D1_R2, D3D2_R2, A_R2, B_R2, C_R2, delta_R2,
                 D1_L1, D2D1_L1, D3D2_L1, A_L1, B_L1, C_L1, delta_L1,
                 D1_R1, D2D1_R1, D3D2_R1, A_R1, B_R1, C_R1, delta_R1,
                 D1_2, D2D1_2, D3D2_2, A_2, B_2, C_2, delta_2, D1_1,
                 D2D1_1, D3D2_1, A_1, B_1, C_1, delta_1)
  
  names(priors) <- c("D1_L2", "D2D1_L2", "D3D2_L2", "A_L2", "B_L2", "C_L2", "delta_L2",
                     "D1_R2", "D2D1_R2", "D3D2_R2", "A_R2", "B_R2", "C_R2", "delta_R2",
                     "D1_L1", "D2D1_L1", "D3D2_L1", "A_L1", "B_L1", "C_L1", "delta_L1",
                     "D1_R1", "D2D1_R1", "D3D2_R1", "A_R1", "B_R1", "C_R1", "delta_R1",
                     "D1_2", "D2D1_2", "D3D2_2", "A_2", "B_2", "C_2", "delta_2", "D1_1",
                     "D2D1_1", "D3D2_1", "A_1", "B_1", "C_1", "delta_1")
  # c("a0_2", "a1_2",
  #   "a2_2", "a3_2", "a4_2", "delta_2", "a0_1", "a1_1", "a2_1", "a3_1",
  #   "a4_1", "delta_1", "CasaleModelIc_2", "CasaleModelIc_1",
  #   "CasaleModelIIa0_2", "CasaleModelIIa1_2", "CasaleModelIIa4_2",
  #   "CasaleModelIIa0_1", "CasaleModelIIa1_1", "CasaleModelIIa4_1",
  #   "CasaleModelIIIa0_2", "CasaleModelIIIa1_2", "CasaleModelIIIa4_2",
  #   "CasaleModelIIIa0_1", "CasaleModelIIIa1_1", "CasaleModelIIIa4_1",
  #   "CasaleModelIVa0_2", "CasaleModelIVa1_2", "CasaleModelIVa2_2",
  #   "CasaleModelIVa3_2", "CasaleModelIVa4_2", "CasaleModelIVa0_1",
  #   "CasaleModelIVa1_1", "CasaleModelIVa2_1", "CasaleModelIVa3_1",
  #   "CasaleModelIVa4_1", "CasaleModelVa0_2", "CasaleModelVa1_2",
  #   "CasaleModelVa2_2", "CasaleModelVa3_2", "CasaleModelVa4_2",
  #   "CasaleModelVa0_1", "CasaleModelVa1_1", "CasaleModelVa2_1",
  #   "CasaleModelVa3_1", "CasaleModelVa4_1")
  
  
  prencours <- NULL
  
  for (i in 1:length(par)) {
    prencours <- c(prencours, priors[[names(par)[i]]])
  }
  
  
  
  parametersMCMC <- matrix(prencours, ncol=7, byrow=T)
  colnames(parametersMCMC) <- c("Density", "Prior1", "Prior2", "SDProp", "Min", "Max", "Init")
  rownames(parametersMCMC)<-names(par)
  parametersMCMC <- as.data.frame(parametersMCMC, stringsAsFactors = FALSE)
  
  for (i in 2:7)
    parametersMCMC[,i] <- as.numeric(parametersMCMC[,i])
  
  parameters <- parametersMCMC
  
  if (accept) {
    return(parameters)
  } else {
    
    repeat {
      
      cat("Proposition:\n")
      print(parameters)
      cat("Name of the parameter to change or Enter to quit:\n")
      f<-scan(nmax=1, quiet=TRUE, what=character())
      
      if (length(f)==0) f <- "q"
      
      if (f=="q") {
        return(parameters)
        
      } else {
        
        variable <- which(f==names(par))
        if (length(variable)==0) {
          cat("The parameter does not exist:\n")
        } else {
          print(variable)
          cat(paste("Change for the parameter ",names(par)[variable],":\n",sep=""))
          
          cat(paste("Distribution of the prior (Enter for default ",parameters[variable, "Density"], "):", sep=""))
          density<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(density)!=0) { parameters[variable, "Density"] <- density } else { density <- parameters[variable, "Density"] }
          
          if (density == "dunif") {
            
            cat(paste("Distribution of the prior, Minimum (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
            f<-scan(nmax=1, quiet=TRUE, what=character())
            if (length(f)!=0) parameters[variable, "Prior1"] <- as.numeric(f)
            cat(paste("Distribution of the prior, Maximum (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
            f<-scan(nmax=1, quiet=TRUE, what=character())
            if (length(f)!=0) parameters[variable, "Prior2"] <- as.numeric(f)
            
          } else {
            
            if (density == "dnorm") {
              
              cat(paste("Distribution of the prior, Mean (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior1"] <- as.numeric(f)
              cat(paste("Distribution of the prior, Standard deviation (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior2"] <- as.numeric(f)
              
            } else {
              
              cat(paste("Distribution of the prior, value 1 (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior1"] <- as.numeric(f)
              cat(paste("Distribution of the prior, value 2 (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior2"] <- as.numeric(f)
              
            }
          }
          
          
          cat(paste("SD of new proposition (Enter for default ",parameters[variable, "SDProp"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "SDProp"] <- as.numeric(f)
          cat(paste("Minimum for the parameter (default ",parameters[variable, "Min"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "Min"] <- as.numeric(f)
          cat(paste("Maximum for the parameter (Enter for default ",parameters[variable, "Max"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "Max"] <- as.numeric(f)
          cat(paste("Initial value (Enter for default ",parameters[variable, "Init"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "Init"] <- as.numeric(f)
        }
        
      }
      
    }
    
  }
  
  for (i in 1:nrow(parameters)) {
    if (parameters[i, "Density"]=="dunif") {
      mn <- max(as.numeric(parameters[i, "Prior1"]), as.numeric(parameters[i, "Min"]))    
    } else {
      mn <- as.numeric(parameters[i, "Min"])
      mx <- as.numeric(parameters[i, "Max"])
    }  
    if (findInterval(as.numeric(parameters[i, "Init"]), c(mn, mx)) != 1) {
      parameters[i, "Init"] <- as.character(mn+(mx-mn)/2)
      warning(paste0("Initial value for parameter ", rownames(parameters)[i], " was out of range; It has been corrected. Check it.")) 
    }
  }
  
  
}
