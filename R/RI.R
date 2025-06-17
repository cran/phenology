#' RI returns an expected remigration interval
#' @title Return an expected remigration interval.
#' @author Marc Girondot
#' @return Return a remigration interval.\cr
#' @param s Time-conditional probability of survival
#' @param t Time-conditional probability of tag retention
#' @param r Time-conditional probability of return
#' @param c Probability of return
#' @param p Annual probability of observation
#' @description Model of remigration interval\cr
#' Note that r, s and t are conditional probabilities. If c is null, then return probabilities are 
#' estimated from r. r can be named vector. For example:\cr
#' r <- c(r1=0.5, r2=0.60, r3=1) is equivalent to c <- c(c1=0.5, c2=0.3, c3=0.2)\cr
#' The vector of r described the probability that a female returned after 
#' 1, 2, 3 years among those who have not nested before. 
#' The vector of r is the same but defining the return probability for an initial female.\cr
#' @family Model of Remigration Interval
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' s <- c(s1=1, s2=1, s3=1, s4=1, s5=1)
#' t <- c(t1=0.95, t2=1, t3=1, t4=1, t5=1) 
#' r <- c(r1=0.1, r2=0.8, r3=0.7, r4=0.7, r5=1)
#' p <- c(p1=0.6, p2=0.6, p3=0.6, p4=0.6, p5=0.6)
#' 
#' # r is equivalent to 
#' c <- c(c1=0.1, c2=0.72, c3=0.126, c4=0.0378, c5=0.0162)
#' # Then the true remigration interval is:
#' ri_true <- sum(1:5*c[1:5])
#' 
#' # BP is then
#' RI2BP(proportion=c, RI=c(mean=1:5, sd=rep(0, 5)))
#' 
#' s_ri  <- NULL
#' for (sx in seq(from=0.01, to=1, by=0.01)) {
#'   s[] <- sx
#'   ri1 <- RI(s=s, t=t, r=r, p=p)
#'   s_ri  <- c(s_ri,sum(1:5*ri1)/sum(ri1))
#' }
#' 
#' par(mar=c(4, 4, 1, 1)+0.4)
#' 
#' plot(x=seq(from=0.01, to=1, by=0.01), y=s_ri, type="l", 
#'      las=1, bty="n", ylim=c(0, 4), 
#'      xlab="Annuual survival probabilities", ylab="Naive Remigration Interval", 
#'     main="")
#' segments(x0=0.01, x1=1, y0=ri_true, y1=ri_true, lty=2, col="red")
#' legend("topright", legend="True remigration interval", lty=2, col="red")
#' 
#' }
#' @export


RI <- function(s, t, r=NULL, c=NULL, p) {
  

  if (names(s[1]) != "s") {
    s <- s[order(as.numeric(substr(names(s), 2, nchar(s))))]
  }
 
  if (names(t[1]) != "t") {
    t <- t[order(as.numeric(substr(names(t), 2, nchar(t))))]
  }
  
  if (!is.null(r)) {
    if (names(r[1]) != "r") {
      r <- r[order(as.numeric(substr(names(r), 2, nchar(r))))]
    }
  }
  
  if (!is.null(c)) {
    if (names(c[1]) != "c") {
      c <- c[order(as.numeric(substr(names(c), 2, nchar(c))))]
    }
  }
  
  if (names(p[1]) != "p") {
    p <- p[order(as.numeric(substr(names(p), 2, nchar(p))))]
  }
  
  
  
  if (is.null(c)) {
    
    rp <- rep(NA, length(r))
    rp[1] <- r[1]
    if (length(r) != 1) {
      for (i in 2:length(r)) rp[i] <- r[i]*prod(1 - r[1:(i-1)])
    }
    r <- rp
  } else {
    r <- c
  }
  
  
  # r <- r/sum(r)
  k <- matrix(data = 1, nrow=1)
  N <- 1

  # Note that r, s and t are conditional probabilities

  n.return.vue <- NULL

  for (l in 1:length(s)) {
    n.survival <- N * prod(s[1:l])
    n.tagretention <- n.survival * prod(t[1:l])
    n.return <- 0
  for (h in 1:nrow(k)) {
    ts <- k[h, , drop=TRUE]
    j <- 1
    l <- 1
    if (length(ts) != 1) {
    for (n in 1:(length(ts)-1)) {
      if (ts[n]==1) {
        j <- j*r[l]*(1-p[n])
        l <- 1
      } else {
        l <- l + 1
      }
    }
    }
    j <- j*r[l]*p[length(s)]
    n.return <- n.return + j * n.tagretention
  }
  n.return.vue <- c(n.return.vue, n.return)
  # je prépare le suivant
k2 <- k
k2[, ncol(k)] <- 0
k3 <- rbind(k, k2)
k4 <- cbind(k3, rep(1, nrow(k3)))
k <- k4
}
return(unname(n.return.vue))
}


