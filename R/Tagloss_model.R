#' Tagloss_model returns the daily rate of tag loss.
#' @title Return the daily rate of tag loss.
#' @author Marc Girondot
#' @return Return the daily rate of tag loss.
#' @param t Time for which values of model must be estimated
#' @param par Parameter
#' @description This function compute a model of daily tag loss rate for days t 
#' based on a set of parameters, par.\cr
#' @family Model of Tag-loss
#' @examples
#' library(phenology)
#' \dontrun{
#' # Example
#' t <- 1:1000
#' par <- c(D1=200, D2D1=100, D3D2=200, 
#'          A=-logit(0.02), B=-logit(0.05), C=-logit(0.07))
#' y <- Tagloss_model(t, par)
#' plot(x=t, y, type="l")
#' }
#' @export

Tagloss_model <- function(t, par) {
  # print(par)
  names(par) <- gsub("_.*", "", names(par))
  if (!is.na(par["D1"])) D1 <- abs(par["D1"]) else D1 <- 0
  if (!is.na(par["D2D1"])) D2D1 <- abs(par["D2D1"]) else D2D1 <- max(t)+1
  D2 <- D1 + D2D1
  if (!is.na(par["D3D2"])) D3D2 <- abs(par["D3D2"]) else D3D2 <- max(t)+1
  D3 <- D2 + D3D2
  if (!is.na(par["A"])) A <- invlogit(-par["A"]) else A <- 0
  if (A == 0) A <- 1E-9
  if (A == 1) A <- 1-1E-9
  if (!is.na(par["B"])) B <- invlogit(-par["B"]) else B <- 0
  if (B == 0) B <- 1E-9
  if (B == 1) B <- 1-1E-9
  if (!is.na(par["C"])) C <- invlogit(-par["C"]) else C <- 0
  if (C == 0) C <- 1E-9
  if (C == 1) C <- 1-1E-9
  
  if (is.na(par["delta"])) delta <- 1 else delta <- invlogit(-par["delta"])
  if (D1==0) D1_2 <- 0.1 else D1_2 <- 2*D1
  if (D3==D2) D3D2 <- (D3-0.1)*D2 else D3D2 <-D3-D2
  
  z <- ifelse(t<=D1, 
         (1+cos(pi*((t+D1)/(D1_2))))*(A-B)+B, 
         ifelse(t<=D2, B, 
                ifelse(t<=D3, (1+cos(pi*((t-D2)/(D3D2))))*0.5*(B-C)+C, C))) * delta
  ifelse(z<0, 1E-9, ifelse(z>1, 1-1E-9, z))
}
