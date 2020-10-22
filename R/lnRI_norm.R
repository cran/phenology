#' lnRI_norm returns a ln L
#' @title Return a remigration interval.
#' @author Marc Girondot
#' @return Return a remigration interval.\cr
#' @param data Data with remigration intervals
#' @param x Vector of parameters
#' @param kl Maximum number of years for remigration intervals.
#' @description Model of remigration interval\cr
#' The vector of parameters must include:\cr
#' sx, survival for year x\cr
#' if s is included, all years have the same survival\cr
#' tx, Tag retention for year x\cr
#' rx, probability of return for year x\cr
#' cx, probability of return for year x\cr
#' px, probability of observation for year x\cr
#' @family Model of Remigration Interval
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' # Each year a fraction of 0.9 is surviving
#' s <- c(s1=0.9, s2=0.9, s3=0.9, s4=0.9, s5=0.9)
#' # Probability of tag retention; 0.95 the first year then after no loss
#' t <- c(t1=0.95, t2=1, t3=1, t4=1, t5=1)
#' # Time-conditional return probability - This is the true remigration rate
#' r <- c(r1=0.1, r2=0.8, r3=0.7, r4=0.7, r5=1)
#' # Capture probability
#' p <- c(p1=0.6, p2=0.6, p3=0.6, p4=0.6, p5=0.6)
#' # Number of observations for 400 tagged females after 1, 2, 3, 4, and 5 years
#' OBS <- c(400, 10, 120, 40, 20, 10)
#' # Likelihood of the observed number based on the model
#' LnRI_norm(data=OBS, x = c(s, t, r, p, sd=2) )
#' LnRI_norm(data=OBS, x = c(s=0.97, t, r, p, sd=2) )
#' 
#' }
#' @export

LnRI_norm <- function(data, x, kl=NULL) {
  
  par_s <- x[substr(names(x), 1, 1)=="s" & substr(names(x), 1, 4)!="sd"]
  if (names(par_s[1]) != "s") {
    par_s <- par_s[order(as.numeric(substr(names(par_s), 2, nchar(par_s))))]
  }
  
  par_t <- x[substr(names(x), 1, 1)=="t"]
  if (names(par_t[1]) != "t") {
    par_t <- par_t[order(as.numeric(substr(names(par_t), 2, nchar(par_t))))]
  }
  
  par_r <- x[substr(names(x), 1, 1)=="r"]
  if (identical(structure(numeric(0), .Names = character(0)), par_r)) {
    par_r <- NULL
  } else {
    if (names(par_r[1]) != "r") {
      par_r <- par_r[order(as.numeric(substr(names(par_r), 2, nchar(par_r))))]
    }
  }
  
  par_c <- x[substr(names(x), 1, 1)=="c"]
  if (identical(structure(numeric(0), .Names = character(0)), par_c)) {
    par_c <- NULL
    } else {
      if (names(par_c[1]) != "c") {
        par_c <- par_c[order(as.numeric(substr(names(par_c), 2, nchar(par_c))))]
      }
    }
  
  par_p <- x[substr(names(x), 1, 1)=="p"]
  if (names(par_p[1]) != "p") {
    par_p <- par_p[order(as.numeric(substr(names(par_p), 2, nchar(par_p))))]
  }
  
  if (is.null(kl)) {
    kl <- max(length(par_s), length(par_t), length(par_c), length(par_r), length(par_p), length(data)-1)
  }
  
  par_s <- c(par_s, rep(par_s[length(par_s)], kl-length(par_s)))
  names(par_s) <- paste0("s", 1:kl)
  par_t <- c(par_t, rep(par_t[length(par_t)], kl-length(par_t)))
  names(par_t) <- paste0("t", 1:kl)
  if (!is.null(par_c)) {
    par_c <- c(par_c, rep(par_c[length(par_c)], kl-length(par_c)))
    names(par_c) <- paste0("c", 1:kl)
  }
  if (!is.null(par_r)) {
    par_r <- c(par_r, rep(par_r[length(par_r)], kl-length(par_r)))
    names(par_r) <- paste0("r", 1:kl)
  }
  
  par_p <- c(par_p, rep(par_p[length(par_p)], kl-length(par_p)))
  names(par_p) <- paste0("p", 1:kl)
  
  sd <- x["sd"]
  n.return.vue <- RI(s=par_s, # survival
                     t=par_t, # Tag loss
                     r=par_r, # probability of return
                     c=par_c, # probability of return
                     p=par_p  # probability of observation
  )
  N0 <- data[1]
  n.return.vue <- n.return.vue*N0
  data <- data[-1]
  
  kl <- min(length(data), kl)
  
  n.return.vue <- n.return.vue[1:kl]
  data <- data[1:kl]
  
  n.return.vue <- c(N0-sum(n.return.vue), n.return.vue)
  data <- c(N0-sum(data), data)
  
  return(-sum(dnorm(data, mean=n.return.vue, 
                    sd = ifelse(n.return.vue*sd<=0, 1, n.return.vue*sd), log=TRUE)))
}

