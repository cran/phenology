# .Lnegbin estimate a negative binomial likelihood.
# @title The function ".Lnegbin"
# @author Marc Girondot
# @return Return the likelihood
# @param x Set of parameters
# @param pt Transfer parameters
# pt=list(data=data, fixed=fixed.parameters, 
#         model_before=model_before, 
#         parallel=parallel, 
#         out=out, 
#         namespar=names(fitted.parameters), 
#         zero=zero, cofactors=cofactors, 
#         add.cofactors=add.cofactors, 
#         WAIC=TRUE, WAIC.bybeach=FALSE)
# @description Function of the package phenology

.Lnegbin <- function (x, pt) {
  
  if (!is.null(pt$parallel)) {
    if (pt$parallel) {
      mc.cores <- getOption("mc.cores", detectCores())
      forking <- getOption("forking", ifelse(.Platform$OS.type == "windows", FALSE, TRUE))
    } else {
      mc.cores <- 1
      forking <- FALSE
    }
  } else {
    mc.cores <- getOption("mc.cores", detectCores())
    forking <- getOption("forking", ifelse(.Platform$OS.type == "windows", FALSE, TRUE))
  }
  
  
  
  # 2/10/2019: A priori ca sert en méthode Brent qui ne transmet pas les noms. Oui c'est ça
  if (is.null(names(x))) {
    if (!is.null(pt$namespar)) {
      names(x) <- pt$namespar
    }
  }
  
  # 19/3/2016: je rajoute cofactors et add.cofactors dans pt
  
  #  if (length(pt)==1 & names(pt)[1]=="pt") pt <- pt$pt
  # .phenology.env<- NULL
  # rm(.phenology.env)
  
  # pt=list(data=data, fixed=fixed.parameters)
  
  sum = 0
  # je mets tous les paramètres dans xpar
  xpar <- c(x, pt$fixed)
  # Si out vaut TRUE, renvoie la somme des vraisemblances
  # Sinon la vraisemblance de chaque série
  
  if (is.null(pt$out)) {
    out <- TRUE
  } else {
    out <- pt$out
  }
  
  
  if (is.null(pt$zero)) {
    zero <- 1E-9
  } else {
    zero <- pt$zero
  }
  
  # deb <- ifelse(zero, 0, 1)
  model_before <- pt$model_before
  
  datatot <- pt$data
  
  add.cofactors <- pt$add.cofactors
  cofactors <- pt$cofactors
  
  
  daily_count <- getFromNamespace(".daily_count", ns = "phenology")
  format_par <- getFromNamespace(".format_par", ns = "phenology")
  
  rg <- try(universalmclapply(X = seq_along(datatot), mc.cores =  mc.cores,
                              FUN = function(k) {
                                
                                # for (k in seq_along(datatot)) {
                                # 
                                # print(k)
                                # k <- 1
                                
                                deb <- 0
                                
                                data <- datatot[[k]]
                                nmser <- names(datatot)[k]
                                
                                xparec <- format_par(xpar=xpar, serie=nmser, model_before=model_before)
                                th <- xparec["Theta"]
                                if (th < 1E-6) th <- 1E-6
                                # Il faut que je fasse daily_count(data$ordinal[i], xparec, print = FALSE, zero = zero) une seule fois
                                # Le premier jour de la série est 0 et donc il est à l'indice 1
                                # data$ordinal[i]
                                # data$ordinal2[i]
                                
                                j <- unname(unlist(apply(X = data, MARGIN = 1, FUN=function(x) seq(from=x["ordinal"], to=ifelse(is.na(x["ordinal2"]), x["ordinal"], x["ordinal2"]), by=1))))
                                nbperday <- rep(NA, max(366, j+1)-deb)
                                # Je dois envoyer les données de cofactors avec seulement les dates
                                # à analyser
                                
                                cof <- NULL
                                if ((!is.null(add.cofactors)) & (!is.null(cofactors))) {
                                  cof <- cofactors[cofactors$Date %in% (data[1, "Date"]+j), ]
                                  cof <- cof[, -1, drop=FALSE]
                                  cof <- as.data.frame(cbind(Date=j, cof))
                                }
                                
                                nbperday[j+1] <- daily_count(d=unname(j), xpar=xparec, 
                                                             cofactors=cof, 
                                                             add.cofactors=add.cofactors, 
                                                             print = FALSE, zero = zero)
                                
                                
                                
                                for (i in 1:nrow(data)) {
                                  
                                  if ((is.na(data[i, "A"])) | (is.na(data[i, "S"])) | (is.null(data[i, "A"])) | (is.null(data[i, "S"]))) {
                                    d <- 1
                                  } else {
                                    d <- 1/(1+exp(-(1/(4*data[i, "S"]))*(data[i, "A"]-0)))
                                  }
                                  
                                  if (data$ZeroCounts[i]) {
                                    
                                    if (data$CountTypes[i] == "exact") {
                                      if (is.na(data$Date2[i])) {
                                        
                                        # sumnbcount <- daily_count(data$ordinal[i], xparec, print = FALSE, zero = zero)
                                        # Le premier jour de la série est 0 et donc il est à l'indice 1
                                        sumnbcount <- nbperday[data$ordinal[i]+1]
                                        
                                        lnli2 <- -dnbinom(data$nombre[i], size = th, 
                                                          mu = sumnbcount*d, log = TRUE)
                                      } else {
                                        # break
                                        nbjour <- data$ordinal2[i] - data$ordinal[i] + 1
                                        # nbcount <- daily_count((1:nbjour) + data$ordinal[i] - 1, xparec, print = FALSE)
                                        # Le premier jour de la série est 0 et donc il est à l'indice 1
                                        # nbcount <- nbperday[((1:nbjour) + data$ordinal[i] - 1)+1]
                                        nbcount <- nbperday[(1:nbjour) + data$ordinal[i]]
                                        sumnbcount <- sum(nbcount)
                                        
                                        if (!is.null(pt$method_Snbinom)) {
                                          method_Snbinom <- pt$method_Snbinom
                                        } else {
                                          method_Snbinom <- "saddlepoint"
                                        }
                                        
                                        if ((is.na(data[i, "A"])) | (is.na(data[i, "S"])) | (is.null(data[i, "A"])) | (is.null(data[i, "S"]))) {
                                          d <- rep(1, nbjour)
                                        } else {
                                          d <- 1/(1+exp(-(1/(4*data[i, "S"]))*(data[i, "A"]-(0:(nbjour-1)))))
                                        }
                                        lnli2 <- -dSnbinom(x=data$nombre[i], 
                                                           size = th, 
                                                           mu = nbcount*rev(d), 
                                                           log = TRUE, 
                                                           method = method_Snbinom, 
                                                           normalize = FALSE)
                                      }
                                    } else {
                                      # je suis en minimum
                                      if (data$CountTypes[i] == "minimum") {
                                        if (is.na(data$Date2[i])) {
                                          
                                          # sumnbcount <- daily_count(data$ordinal[i], xparec, print = FALSE, zero = zero)
                                          # Le premier jour de la série est 0 et donc il est à l'indice 1
                                          sumnbcount <- nbperday[data$ordinal[i]+1]
                                          # Par exemple j'ai un minimum de 1
                                          # pnbinom(q=1-1, size = 3, mu = 10, log.p = FALSE, lower.tail = FALSE)
                                          # C'est bon puisque 
                                          # pnbinom(q=1-1, size = 3, mu = 10, log.p = FALSE, lower.tail = FALSE)+dnbinom(0, size = 3, mu = 10, log=FALSE) = 1
                                          if (data$nombre[i]-1 >= 0) { 
                                            
                                            lnli2 <- -pnbinom(q=data$nombre[i]-1, size = th, 
                                                              mu = sumnbcount*d, log.p = TRUE, 
                                                              lower.tail = TRUE)
                                            
                                            
                                          } else {
                                            # si j'ai 0 observation forcément likelihood=1
                                            lnli2 <- 0
                                          }
                                        } else {
                                          # j'ai une range de dates
                                          nbjour <- data$ordinal2[i] - data$ordinal[i] + 1
                                          # nbcount <- daily_count((1:nbjour) + data$ordinal[i] - 1, xparec, print = FALSE)
                                          # Le premier jour de la série est 0 et donc il est à l'indice 1
                                          nbcount <- nbperday[((1:nbjour) + data$ordinal[i] - 1)+1]
                                          sumnbcount <- sum(nbcount)
                                          if (data$nombre[i]-1 >= 0) { 
                                            if (!is.null(pt$method_Snbinom)) {
                                              method_Snbinom <- pt$method_Snbinom
                                            } else {
                                              method_Snbinom <- "saddlepoint"
                                            }
                                            # -log((1-sum(dSnbinom(x=0:(data$nombre[i]-1), size = th, mu = nbcount, log=FALSE, method = method_Snbinom))))
                                            if ((is.na(data[i, "A"])) | (is.na(data[i, "S"])) | (is.null(data[i, "A"])) | (is.null(data[i, "S"]))) {
                                              d <- rep(1, nbjour)
                                            } else {
                                              d <- 1/(1+exp(-(1/(4*data[i, "S"]))*(data[i, "A"]-(0:(nbjour-1)))))
                                            }
                                            lnli2 <- -pSnbinom(q=data$nombre[i]-1, size = th, mu = nbcount*rev(d), log.p = TRUE, 
                                                               lower.tail = FALSE, method = method_Snbinom, tol=1E-6, normalize = FALSE)
                                            
                                            
                                          } else {
                                            lnli2 <- 0
                                          }
                                        }
                                      } else {
                                        # Je suis en nombre max
                                        nbmin <- data$nombre[i]
                                        nbmax <- as.numeric(data$CountTypes[i])
                                        
                                        if (is.na(data$Date2[i])) {
                                          
                                          # sumnbcount <- daily_count(data$ordinal[i], xparec, print = FALSE, zero = zero)
                                          # Le premier jour de la série est 0 et donc il est à l'indice 1
                                          sumnbcount <- nbperday[data$ordinal[i]+1]
                                          # Par exemple j'ai un minimum de 1
                                          # pnbinom(q=1-1, size = 3, mu = 10, log.p = FALSE, lower.tail = FALSE)
                                          # C'est bon puisque 
                                          # pnbinom(q=1-1, size = 3, mu = 10, log.p = FALSE, lower.tail = FALSE)+dnbinom(0, size = 3, mu = 10, log=FALSE) = 1
                                          
                                          lnli2 <- -sum(dnbinom(x=nbmin:nbmax, size = th, 
                                                                mu = sumnbcount*d, log = TRUE))
                                          
                                        } else {
                                          # j'ai une range de dates
                                          nbjour <- data$ordinal2[i] - data$ordinal[i] + 1
                                          # nbcount <- daily_count((1:nbjour) + data$ordinal[i] - 1, xparec, print = FALSE)
                                          # Le premier jour de la série est 0 et donc il est à l'indice 1
                                          nbcount <- nbperday[((1:nbjour) + data$ordinal[i] - 1)+1]
                                          sumnbcount <- sum(nbcount)
                                          if (data$nombre[i]-1 >= 0) { 
                                            if (!is.null(pt$method_Snbinom)) {
                                              method_Snbinom <- pt$method_Snbinom
                                            } else {
                                              method_Snbinom <- "saddlepoint"
                                            }
                                            # -log((1-sum(dSnbinom(x=0:(data$nombre[i]-1), size = th, mu = nbcount, log=FALSE, method = method_Snbinom))))
                                            if ((is.na(data[i, "A"])) | (is.na(data[i, "S"])) | (is.null(data[i, "A"])) | (is.null(data[i, "S"]))) {
                                              d <- rep(1, nbjour)
                                            } else {
                                              d <- 1/(1+exp(-(1/(4*data[i, "S"]))*(data[i, "A"]-(0:(nbjour-1)))))
                                            }
                                            lnli2 <- -sum(dSnbinom(x=nbmin:nbmax, size = th, mu = nbcount*rev(d), log = TRUE, 
                                                                   method = method_Snbinom, tol=1E-6, normalize = FALSE))
                                            
                                          } else {
                                            lnli2 <- 0
                                          }
                                        }
                                      }
                                    }
                                    
                                  } else {
                                    # Là je suis en ZeroCounts == FALSE
                                    
                                    if (data$CountTypes[i] == "exact") {
                                      if (is.na(data$Date2[i])) {
                                        
                                        # sumnbcount <- daily_count(data$ordinal[i], xparec, print = FALSE, zero = zero)
                                        # Le premier jour de la série est 0 et donc il est à l'indice 1
                                        sumnbcount <- nbperday[data$ordinal[i]+1]
                                        
                                        
                                        lnli2 <- (- dnbinom(data$nombre[i], size = th, mu = sumnbcount*d, log = TRUE)+
                                                    log(1 - dnbinom(0, size = th, mu = sumnbcount*d,  log = FALSE)))
                                        
                                      } else {
                                        nbjour <- data$ordinal2[i] - data$ordinal[i] + 1
                                        # nbcount <- daily_count((1:nbjour) + data$ordinal[i] - 1, xparec, print = FALSE)
                                        # Le premier jour de la série est 0 et donc il est à l'indice 1
                                        nbcount <- nbperday[((1:nbjour) + data$ordinal[i] - 1)+1]
                                        sumnbcount <- sum(nbcount)
                                        
                                        # Est ce que les 0 ont été comptés
                                        if (!is.null(pt$method_Snbinom)) {
                                          method_Snbinom <- pt$method_Snbinom
                                        } else {
                                          method_Snbinom <- "furman"
                                          if ((th > 10) | (any(nbcount <1E-5))) {
                                            if (length(nbcount) <= 5) {
                                              method_Snbinom <- "Vellaisamy&Upadhye"
                                            } else {
                                              method_Snbinom <- "approximate.negativebinomial"
                                            }
                                          }
                                          if ((length(nbcount) > 5) & (method_Snbinom == "furman")) {
                                            method_Snbinom <- "saddlepoint"
                                          }
                                        }
                                        if ((is.na(data[i, "A"])) | (is.na(data[i, "S"])) | (is.null(data[i, "A"])) | (is.null(data[i, "S"]))) {
                                          d <- rep(1, nbjour)
                                        } else {
                                          d <- 1/(1+exp(-(1/(4*data[i, "S"]))*(data[i, "A"]-(0:(nbjour-1)))))
                                        }
                                        lnli2 <- (- dSnbinom(data$nombre[i], size = th, mu = nbcount*rev(d), log = TRUE, method = method_Snbinom, normalize = FALSE) + 
                                                    log(1 - dSnbinom(0, size = th, mu = nbcount*rev(d), log = FALSE, method = method_Snbinom, normalize = FALSE)))
                                        
                                      }
                                    } else {
                                      # Je suis en minimum
                                      if (data$CountTypes[i] == "minimum") {
                                        if (is.na(data$Date2[i])) {
                                          
                                          
                                          # sumnbcount <- daily_count(data$ordinal[i], xparec, print = FALSE, zero = zero)
                                          # Le premier jour de la série est 0 et donc il est à l'indice 1
                                          sumnbcount <- nbperday[data$ordinal[i]+1]
                                          
                                          if (data$nombre[i]-1 >= 0) {
                                            
                                            prob.sup.egal.obs <- pnbinom(q=data$nombre[i]-1, size = th, mu = sumnbcount*d, log.p = TRUE, lower.tail = FALSE)
                                            prob0 <- dnbinom(0, size = th, mu = sumnbcount*d,  log = FALSE)
                                            
                                            lnli2 <- (- prob.sup.egal.obs + log(1 - prob0))
                                          } else {
                                            lnli2 <- 0
                                          }
                                        } else {
                                          nbjour <- data$ordinal2[i] - data$ordinal[i] + 1
                                          # nbcount <- daily_count((1:nbjour) + data$ordinal[i] - 1, xparec, print = FALSE)
                                          # Le premier jour de la série est 0 et donc il est à l'indice 1
                                          nbcount <- nbperday[((1:nbjour) + data$ordinal[i] - 1)+1]
                                          sumnbcount <- sum(nbcount)
                                          
                                          # Est ce que les 0 ont été comptés
                                          if (data$nombre[i]-1 >= 0) {
                                            if (!is.null(pt$method_Snbinom)) {
                                              method_Snbinom <- pt$method_Snbinom
                                            } else {
                                              method_Snbinom <- "saddlepoint"
                                            }
                                            if ((is.na(data[i, "A"])) | (is.na(data[i, "S"])) | (is.null(data[i, "A"])) | (is.null(data[i, "S"]))) {
                                              d <- rep(1, nbjour)
                                            } else {
                                              d <- 1/(1+exp(-(1/(4*data[i, "S"]))*(data[i, "A"]-(0:(nbjour-1)))))
                                            }
                                            prob.sup.egal.obs <- pSnbinom(q=data$nombre[i]-1, size = th, mu = nbcount*rev(d), log.p = TRUE, 
                                                                          lower.tail = FALSE, method = method_Snbinom, tol=1E-6, normalize = FALSE)
                                            prob0 <- dSnbinom(0, size = th, mu = nbcount*rev(d),  log = FALSE, method = method_Snbinom, normalize = FALSE)
                                            
                                            lnli2 <- (- prob.sup.egal.obs + log(1 - prob0))
                                            
                                          } else {
                                            lnli2 <- 0
                                          }
                                          
                                        }
                                      } else {
                                        # Je suis en nombre max
                                        nbmin <- data$nombre[i]
                                        nbmax <- as.numeric(data$CountTypes[i])
                                        
                                        if (is.na(data$Date2[i])) {
                                          
                                          
                                          # sumnbcount <- daily_count(data$ordinal[i], xparec, print = FALSE, zero = zero)
                                          # Le premier jour de la série est 0 et donc il est à l'indice 1
                                          sumnbcount <- nbperday[data$ordinal[i]+1]
                                          
                                          if (data$nombre[i]-1 >= 0) {
                                            
                                            prob.sup.egal.obs <- -sum(dnbinom(x=nbmin:nbmax, size = th, 
                                                                              mu = sumnbcount*d, log = TRUE))
                                            
                                            prob0 <- dnbinom(0, size = th, mu = sumnbcount*d,  log = FALSE)
                                            
                                            lnli2 <- (prob.sup.egal.obs + log(1 - prob0))
                                          } else {
                                            lnli2 <- 0
                                          }
                                        } else {
                                          nbjour <- data$ordinal2[i] - data$ordinal[i] + 1
                                          # nbcount <- daily_count((1:nbjour) + data$ordinal[i] - 1, xparec, print = FALSE)
                                          # Le premier jour de la série est 0 et donc il est à l'indice 1
                                          nbcount <- nbperday[((1:nbjour) + data$ordinal[i] - 1)+1]
                                          sumnbcount <- sum(nbcount)
                                          
                                          # Est ce que les 0 ont été comptés
                                          if (data$nombre[i]-1 >= 0) {
                                            if (!is.null(pt$method_Snbinom)) {
                                              method_Snbinom <- pt$method_Snbinom
                                            } else {
                                              method_Snbinom <- "saddlepoint"
                                            }
                                            if ((is.na(data[i, "A"])) | (is.na(data[i, "S"])) | (is.null(data[i, "A"])) | (is.null(data[i, "S"]))) {
                                              d <- rep(1, nbjour)
                                            } else {
                                              d <- 1/(1+exp(-(1/(4*data[i, "S"]))*(data[i, "A"]-(0:(nbjour-1)))))
                                            }
                                            prob.sup.egal.obs <- -sum(dSnbinom(x=nbmin:nbmax, size = th, mu = nbcount*rev(d), log = TRUE, 
                                                                               method = method_Snbinom, tol=1E-6, normalize = FALSE))
                                            prob0 <- dSnbinom(0, size = th, mu = nbcount*rev(d),  log = FALSE, method = method_Snbinom, normalize = FALSE)
                                            
                                            lnli2 <- (prob.sup.egal.obs + log(1 - prob0))
                                            
                                          } else {
                                            lnli2 <- 0
                                          }
                                          
                                        }
                                      }
                                      
                                      
                                    }
                                  }
                                  # print(paste0(as.character(i), "  ", lnli2))
                                  
                                  if (is.na(lnli2)) lnli2 <- 1E9
                                  
                                  data$LnL[i] <- lnli2
                                  data$Modeled[i] <- sumnbcount
                                  
                                }
                                
                                data$LnL[is.infinite(data$LnL)] <- 115.1293 # -log(1e-50)
                                
                                
                                
                                return(data)
                                
                              },
                              clusterExport = list(varlist=c("xpar", "datatot",
                                                             "daily_count", "format_par",
                                                             "zero",
                                                             "model_before", 
                                                             "add.cofactors", "cofactors"),
                                                   envir=environment()),
                              clusterEvalQ=list(expr=expression(library(phenology))
                              ),
                              forking = forking), silent = FALSE)
  
  # }
  
  if (inherits(rg, "try-error")) {
    save(x, file = "x.Rdata")
    save(pt, file = "pt.Rdata")
    stop("Error during likelihood estimation; look at x.Rdata and pt.Rdata")
  }
  
  sum <- sum(sapply(X = rg, function(x) sum(x$LnL, na.rm = TRUE)))
  
  if (!is.null(pt$WAIC)) {
    if (pt$WAIC) {
      LnLi <- lapply(X = rg, function(x) x$LnL)
      if (!is.null(pt$WAIC.bybeach)) {
        if (pt$WAIC.bybeach) {
          LnLi <- unlist(lapply(X = LnLi, function(x) sum(x, na.rm = TRUE)))
        } else {
          LnLi <- unlist(LnLi)
        }
      } else {
        LnLi <- unlist(LnLi)
      }
      
      attributes(sum) <- list(WAIC=-LnLi)
    }
  }
  
  if (is.na(sum)) sum <- 1E9
  
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
