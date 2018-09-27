# .Lnegbin estimate a negative binomial likelihood.
# @title The function ".Lnegbin"
# @author Marc Girondot
# @return Return the likelihood
# @param x Set of parameters
# @param pt Transfer parameters
# @description Function of the package phenology

.Lnegbin <- function (x, pt) {
  
    windows <- (.Platform$OS.type == "windows")
    
    # pt$parallel = FALSE if no parallel (but mclapply is still used)
    mcutil <- 1
    if (!is.null(pt$parallel)) {
      if (pt$parallel) {
        mcutil <- detectCores()
      } 
    } else {
      mcutil <- detectCores()
    }
    
    if (!is.null(pt$namespar)) {
      names(x) <- pt$namespar
    }

  
  # 19/3/2016: je rajoute cofactors et add.cofactors dans pt
  
  #  if (length(pt)==1 & names(pt)[1]=="pt") pt <- pt$pt
  # .phenology.env<- NULL
  # rm(.phenology.env)
  
  # pt=list(data=data, fixed=fixed.parameters, incertitude=method_incertitude, zerocounts=zero_counts)
  
    sum = 0
  # je mets tous les paramètres dans xpar
    xpar <- c(x, pt$fixed)
  # Si out vaut TRUE, renvoie la somme des vraisemblances
  # Sinon la vraisemblance de chaque série
    out <- pt$out
    
    datatot <- pt$data
    infinite <- pt$infinite
    daily_count <- getFromNamespace(".daily_count", ns = "phenology")
    format_par <- getFromNamespace(".format_par", ns = "phenology")
    
    if (mcutil != 1) {
    
    
    if (windows) {

        cl <- makeCluster(mcutil)
        clusterEvalQ(cl = cl, library(phenology))
      
        clusterExport(cl = cl, varlist = c("datatot", "daily_count", 
            "format_par", "pt", "xpar"), envir = environment())
        rg <- parLapplyLB(cl = cl, X = seq_along(datatot), fun = function(k) {
            data <- datatot[[k]]
            nmser <- names(datatot)[k]
            zero_counts <- pt$zerocounts[k]
            zero <- pt$zero
            deb <- ifelse(zero, 0, 1)
            xparec <- format_par(xpar, nmser)
            th <- xparec["Theta"]
            for (i in 1:nrow(data)) {
                if ((data$nombre[i] != 0) || zero_counts) {
                  if (is.na(data$Date2[i])) {
                    sumnbcount <- daily_count(data$ordinal[i], 
                      xparec, print = FALSE, zero = pt$zero)
                    if (!is.null(pt$add.cofactors)) {
                      sumnbcount <- sumnbcount + sum(pt$cofactors[pt$cofactors$Date == 
                        data$Date[i], pt$add.cofactor] * xparec[pt$add.cofactor])
                      if (sumnbcount <= pt$zero) 
                        sumnbcount <- pt$zero
                    }
                    if (!zero_counts) {
                      lnli2 <- -log(dnbinom(data$nombre[i], size = th, 
                        mu = sumnbcount, log = FALSE)/(1 - dnbinom(0, 
                        size = th, mu = sumnbcount, log = FALSE)))
                    } else {
                      lnli2 <- -dnbinom(data$nombre[i], size = th, 
                        mu = sumnbcount, log = TRUE)
                    }
                  } else {
                    nbjour <- data$ordinal2[i] - data$ordinal[i] + 1
                    nbcount <- daily_count((1:nbjour) + data$ordinal[i] - 1, xparec, print = FALSE)
                    sumnbcount <- sum(nbcount)
                    if (pt$incertitude == 1 | pt$incertitude == 
                      "convolution") {
                      if (!zero_counts) {
                        lnli2 <- -log(dSnbinom(data$nombre[i], 
                          size = th, mu = nbcount, log = FALSE, 
                          infinite = infinite)/(1 - dSnbinom(0, 
                          size = th, mu = nbcount, log = FALSE, 
                          infinite = infinite)))
                      }                      else {
                        lnli2 <- -dSnbinom(data$nombre[i], size = th, 
                          mu = nbcount, log = TRUE, infinite = infinite)
                      }
                    }
                    if (pt$incertitude == 2 | pt$incertitude == "combinatory") {
                      nbcountrel <- nbcount/sumnbcount
                      nbcountrel[nbcountrel == 0] <- pt$zero
                      nbcountrel[nbcountrel == 1] <- 1 - pt$zero
                      N <- data$nombre[i]
                      xx <- combn(N + nbjour - 1, nbjour - 1)
                      a <- cbind(0, diag(nbjour)) - cbind(diag(nbjour), 
                        0)
                      tb <- t(a %*% rbind(0, xx, N + nbjour) - 
                        1)
                      a <- try(matrix(rep(0, (N + 1) * nbjour), 
                        nrow = N + 1, ncol = nbjour), silent = TRUE)
                      if (class(a) == "try-error") {
                        stop("Too many incertitudes on the days. Use the other method named 'convolution'.")
                      }
                      for (ii in deb:N) {
                        for (countday in 1:nbjour) {
                          if (!zero_counts) {
                            a[ii + 1, countday] <- log(dnbinom(ii, 
                              size = th, mu = nbcount[countday], 
                              log = FALSE)/(1 - dnbinom(0, size = th, 
                              mu = nbcount[countday], log = FALSE)))
                          }                          else {
                            a[ii + 1, countday] <- dnbinom(ii, 
                              size = th, mu = nbcount[countday], 
                              log = TRUE)
                          }
                        }
                      }
                      sump <- 0
                      for (ii in 1:dim(tb)[1]) {
                        p <- dmultinom(tb[ii, 1:nbjour], prob = nbcountrel, 
                          log = TRUE)
                        for (countday in 1:nbjour) {
                          p <- p + a[tb[ii, countday] + 1, countday]
                        }
                        sump <- sump + exp(p)
                      }
                      lnli2 <- -log(sump)
                    }
                  }
                }                else {
                  lnli2 <- NA
                  sumnbcount <- NA
                }
                datatot[[k]]$LnL[i] <- lnli2
                datatot[[k]]$Modeled[i] <- sumnbcount
            }
            return(datatot[[k]])
        })
        stopCluster(cl)
    }    else {
        rg <- mclapply(X = seq_along(datatot), mc.cores =  mcutil, 
            FUN = function(k) {
                data <- datatot[[k]]
                nmser <- names(datatot)[k]
                zero_counts <- pt$zerocounts[k]
                zero <- pt$zero
                deb <- ifelse(zero, 0, 1)
                xparec <- format_par(xpar, nmser)
                th <- xparec["Theta"]
                for (i in 1:nrow(data)) {
                  if ((data$nombre[i] != 0) || zero_counts) {
                    if (is.na(data$Date2[i])) {
                      sumnbcount <- daily_count(data$ordinal[i], 
                        xparec, print = FALSE, zero = pt$zero)
                      if (!is.null(pt$add.cofactors)) {
                        sumnbcount <- sumnbcount + sum(pt$cofactors[pt$cofactors$Date == 
                          data$Date[i], pt$add.cofactor] * xparec[pt$add.cofactor])
                        if (sumnbcount <= pt$zero) 
                          sumnbcount <- pt$zero
                      }
                      if (!zero_counts) {
                        lnli2 <- -log(dnbinom(data$nombre[i], 
                          size = th, mu = sumnbcount, log = FALSE)/(1 - 
                          dnbinom(0, size = th, mu = sumnbcount, 
                            log = FALSE)))
                      }                      else {
                        lnli2 <- -dnbinom(data$nombre[i], size = th, 
                          mu = sumnbcount, log = TRUE)
                      }
                    }                    else {
                      nbjour <- data$ordinal2[i] - data$ordinal[i] + 
                        1
                      nbcount <- daily_count((1:nbjour) + data$ordinal[i] - 
                        1, xparec, print = FALSE)
                      sumnbcount <- sum(nbcount)
                      if (pt$incertitude == 1 | pt$incertitude == 
                        "convolution") {
                        if (!zero_counts) {
                          lnli2 <- -log(dSnbinom(data$nombre[i], 
                            size = th, mu = nbcount, log = FALSE, 
                            infinite = infinite)/(1 - dSnbinom(0, 
                            size = th, mu = nbcount, log = FALSE, 
                            infinite = infinite)))
                        }                        else {
                          lnli2 <- -dSnbinom(data$nombre[i], 
                            size = th, mu = nbcount, log = TRUE, 
                            infinite = infinite)
                        }
                      }
                      if (pt$incertitude == 2 | pt$incertitude == 
                        "combinatory") {
                        nbcountrel <- nbcount/sumnbcount
                        nbcountrel[nbcountrel == 0] <- pt$zero
                        nbcountrel[nbcountrel == 1] <- 1 - pt$zero
                        N <- data$nombre[i]
                        xx <- combn(N + nbjour - 1, nbjour - 
                          1)
                        a <- cbind(0, diag(nbjour)) - cbind(diag(nbjour), 
                          0)
                        tb <- t(a %*% rbind(0, xx, N + nbjour) - 
                          1)
                        a <- try(matrix(rep(0, (N + 1) * nbjour), 
                          nrow = N + 1, ncol = nbjour), silent = TRUE)
                        if (class(a) == "try-error") {
                          stop("Too many incertitudes on the days. Use the other method named 'convolution'.")
                        }
                        for (ii in deb:N) {
                          for (countday in 1:nbjour) {
                            if (!zero_counts) {
                              a[ii + 1, countday] <- log(dnbinom(ii, 
                                size = th, mu = nbcount[countday], 
                                log = FALSE)/(1 - dnbinom(0, 
                                size = th, mu = nbcount[countday], 
                                log = FALSE)))
                            }                            else {
                              a[ii + 1, countday] <- dnbinom(ii, 
                                size = th, mu = nbcount[countday], 
                                log = TRUE)
                            }
                          }
                        }
                        sump <- 0
                        for (ii in 1:dim(tb)[1]) {
                          p <- dmultinom(tb[ii, 1:nbjour], prob = nbcountrel, 
                            log = TRUE)
                          for (countday in 1:nbjour) {
                            p <- p + a[tb[ii, countday] + 1, 
                              countday]
                          }
                          sump <- sump + exp(p)
                        }
                        lnli2 <- -log(sump)
                      }
                    }
                  }                  else {
                    lnli2 <- NA
                    sumnbcount <- NA
                  }
                  datatot[[k]]$LnL[i] <- lnli2
                  datatot[[k]]$Modeled[i] <- sumnbcount
                }
                datatot[[k]]$LnL[is.infinite(datatot[[k]]$LnL)] <- -log(1e-50)
                return(datatot[[k]])
            })
    }
    } else {
      # version non parallele
      
      rg <- lapply(X = seq_along(datatot), 
                     FUN = function(k) {
                       data <- datatot[[k]]
                       nmser <- names(datatot)[k]
                       zero_counts <- pt$zerocounts[k]
                       zero <- pt$zero
                       deb <- ifelse(zero, 0, 1)
                       xparec <- format_par(xpar, nmser)
                       th <- xparec["Theta"]
                       for (i in 1:nrow(data)) {
                         if ((data$nombre[i] != 0) || zero_counts) {
                           if (is.na(data$Date2[i])) {
                             sumnbcount <- daily_count(data$ordinal[i], 
                                                       xparec, print = FALSE, zero = pt$zero)
                             if (!is.null(pt$add.cofactors)) {
                               sumnbcount <- sumnbcount + sum(pt$cofactors[pt$cofactors$Date == 
                                                                             data$Date[i], pt$add.cofactor] * xparec[pt$add.cofactor])
                               if (sumnbcount <= pt$zero) 
                                 sumnbcount <- pt$zero
                             }
                             if (!zero_counts) {
                               lnli2 <- -log(dnbinom(data$nombre[i], 
                                                     size = th, mu = sumnbcount, log = FALSE)/(1 - 
                                                                                                 dnbinom(0, size = th, mu = sumnbcount, 
                                                                                                         log = FALSE)))
                             }                      else {
                               lnli2 <- -dnbinom(data$nombre[i], size = th, 
                                                 mu = sumnbcount, log = TRUE)
                             }
                           }                    else {
                             nbjour <- data$ordinal2[i] - data$ordinal[i] + 
                               1
                             nbcount <- daily_count((1:nbjour) + data$ordinal[i] - 
                                                      1, xparec, print = FALSE)
                             sumnbcount <- sum(nbcount)
                             if (pt$incertitude == 1 | pt$incertitude == 
                                 "convolution") {
                               if (!zero_counts) {
                                 lnli2 <- -log(dSnbinom(data$nombre[i], 
                                                        size = th, mu = nbcount, log = FALSE, 
                                                        infinite = infinite)/(1 - dSnbinom(0, 
                                                                                           size = th, mu = nbcount, log = FALSE, 
                                                                                           infinite = infinite)))
                               }                        else {
                                 lnli2 <- -dSnbinom(data$nombre[i], 
                                                    size = th, mu = nbcount, log = TRUE, 
                                                    infinite = infinite)
                               }
                             }
                             if (pt$incertitude == 2 | pt$incertitude == 
                                 "combinatory") {
                               nbcountrel <- nbcount/sumnbcount
                               nbcountrel[nbcountrel == 0] <- pt$zero
                               nbcountrel[nbcountrel == 1] <- 1 - pt$zero
                               N <- data$nombre[i]
                               xx <- combn(N + nbjour - 1, nbjour - 
                                             1)
                               a <- cbind(0, diag(nbjour)) - cbind(diag(nbjour), 
                                                                   0)
                               tb <- t(a %*% rbind(0, xx, N + nbjour) - 
                                         1)
                               a <- try(matrix(rep(0, (N + 1) * nbjour), 
                                               nrow = N + 1, ncol = nbjour), silent = TRUE)
                               if (class(a) == "try-error") {
                                 stop("Too many incertitudes on the days. Use the other method named 'convolution'.")
                               }
                               for (ii in deb:N) {
                                 for (countday in 1:nbjour) {
                                   if (!zero_counts) {
                                     a[ii + 1, countday] <- log(dnbinom(ii, 
                                                                        size = th, mu = nbcount[countday], 
                                                                        log = FALSE)/(1 - dnbinom(0, 
                                                                                                  size = th, mu = nbcount[countday], 
                                                                                                  log = FALSE)))
                                   }                            else {
                                     a[ii + 1, countday] <- dnbinom(ii, 
                                                                    size = th, mu = nbcount[countday], 
                                                                    log = TRUE)
                                   }
                                 }
                               }
                               sump <- 0
                               for (ii in 1:dim(tb)[1]) {
                                 p <- dmultinom(tb[ii, 1:nbjour], prob = nbcountrel, 
                                                log = TRUE)
                                 for (countday in 1:nbjour) {
                                   p <- p + a[tb[ii, countday] + 1, 
                                              countday]
                                 }
                                 sump <- sump + exp(p)
                               }
                               lnli2 <- -log(sump)
                             }
                           }
                         }                  else {
                           lnli2 <- NA
                           sumnbcount <- NA
                         }
                         datatot[[k]]$LnL[i] <- lnli2
                         datatot[[k]]$Modeled[i] <- sumnbcount
                       }
                       datatot[[k]]$LnL[is.infinite(datatot[[k]]$LnL)] <- -log(1e-50)
                       return(datatot[[k]])
                     })
      
      
    }
      
      
    
      
    sum <- sum(sapply(X = rg, function(x) sum(x$LnL, na.rm = TRUE)))
    if (is.infinite(sum)) {
        save(x, file = "x.Rdata")
        save(pt, file = "pt.Rdata")
    }
    
    if (!is.null(pt$store.intermediate)) {
    if (pt$store.intermediate) {
      load(file=pt$file.intermediate)
      store <- c(store, iteration=list(list(fitted=x, fixed=pt$fixed, loglik=sum)))
      save(store, file=pt$file.intermediate)
    }
}
      
    
    if (out) {
        return(sum)
    }    else {
        names(rg) <- names(datatot)
        return(rg)
    }
}
