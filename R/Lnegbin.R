# .Lnegbin estimate a negative binomial likelihood.
# @title The function ".Lnegbin"
# @author Marc Girondot
# @return Return the likelihood
# @param x Set of parameters
# @param pt Transfer parameters
# pt=list(data=data, fixed=fixed.parameters, 
#         parallel=parallel, 
#         zerocounts=zero_counts, 
#         tol=tol, out=out, 
#         namespar=names(fitted.parameters), 
#         zero=zero, cofactors=cofactors, 
#         add.cofactors=add.cofactors)
# @description Function of the package phenology

.Lnegbin <- function (x, pt) {
  
  # pt$parallel = FALSE if no parallel (but mclapply is still used)
  mcutil <- 1
  if (!is.null(pt$parallel)) {
    if (pt$parallel) {
      mcutil <- detectCores()
    } 
  } else {
    mcutil <- detectCores()
  }
  
  # 12/10/2018. J'ai eu une erreur Error in if (mcutil != 1) { : valeur manquante là où TRUE / FALSE est requis
  if (is.null(mcutil)) mcutil <- 1
  if (is.na(mcutil)) mcutil <- 1
  
  # 2/10/2019: A priori ca sert en méthode Brent qui ne transmet pas les noms. Oui c'est ça
  if (!is.null(pt$namespar)) {
    names(x) <- pt$namespar
  }
  
  
  # 19/3/2016: je rajoute cofactors et add.cofactors dans pt
  
  #  if (length(pt)==1 & names(pt)[1]=="pt") pt <- pt$pt
  # .phenology.env<- NULL
  # rm(.phenology.env)
  
  # pt=list(data=data, fixed=fixed.parameters, zerocounts=zero_counts)
  
  sum = 0
  # je mets tous les paramètres dans xpar
  xpar <- c(x, pt$fixed)
  # Si out vaut TRUE, renvoie la somme des vraisemblances
  # Sinon la vraisemblance de chaque série
  out <- pt$out
  
  datatot <- pt$data
  tol <- pt$tol
  if (is.null(tol)) tol <- 1e-6
  daily_count <- getFromNamespace(".daily_count", ns = "phenology")
  format_par <- getFromNamespace(".format_par", ns = "phenology")
  
  rg <- try(universalmclapply(X = seq_along(datatot), mc.cores =  mcutil, 
                              FUN = function(k) {
                                
                                # for (k in seq_along(datatot)) {
                                #   
                                #   print(k)
                                
                                data <- datatot[[k]]
                                nmser <- names(datatot)[k]
                                zero_counts <- pt$zerocounts[k]
                                zero <- pt$zero
                                deb <- ifelse(zero, 0, 1)
                                xparec <- format_par(xpar, nmser)
                                th <- xparec["Theta"]
                                
                                # Il faut que je fasse daily_count(data$ordinal[i], xparec, print = FALSE, zero = pt$zero) une seule fois
                                # Le premier jour de la série est 0 et donc il est à l'indice 1
                                # data$ordinal[i]
                                # data$ordinal2[i]
                                
                                j <- unlist(apply(X = data, MARGIN = 1, FUN=function(x) x["ordinal"]:ifelse(is.na(x["ordinal2"]), x["ordinal"], x["ordinal2"])))
                                nbperday <- rep(NA, max(366, j+1))
                                # Je dois envoyer les données de cofactors avec seulement les dates
                                # à analyser
                                cof <- NULL
                                if ((!is.null(pt$add.cofactors)) & (!is.null(pt$cofactors))) {
                                  cof <- pt$cofactors[pt$cofactors$Date %in% (data[1, "Date"]+j), ]
                                  cof <- cof[, -1, drop=FALSE]
                                  cof <- as.data.frame(cbind(Date=j, cof))
                                }
                                
                                nbperday[j+1] <- daily_count(d=j, xpar=xparec, 
                                                             cofactors=cof, 
                                                             add.cofactors=pt$add.cofactors, 
                                                             print = FALSE, zero = pt$zero)
                                
                                
                                
                                
                                # if (!is.null(pt$add.cofactors)) {
                                #   # J'ai des cofacteurs
                                #   
                                #   if (zero_counts) {
                                #     
                                #     for (i in 1:nrow(data)) {
                                #       
                                #       if (is.na(data$Date2[i])) {
                                #         
                                #         # sumnbcount <- daily_count(data$ordinal[i], xparec, print = FALSE, zero = pt$zero)
                                #         # Le premier jour de la série est 0 et donc il est à l'indice 1
                                #         sumnbcount <- nbperday[data$ordinal[i]+1]
                                #         
                                #         
                                #         sumnbcount <- sumnbcount + sum(pt$cofactors[pt$cofactors$Date == 
                                #                                                       data$Date[i], pt$add.cofactor] * 
                                #                                          xparec[pt$add.cofactor])
                                #         if (sumnbcount <= pt$zero) sumnbcount <- pt$zero
                                #         # dni <- data$nombre[i]
                                #         lnli2 <- -dnbinom(data$nombre[i], size = th, 
                                #                            mu = sumnbcount, log = TRUE)
                                #         
                                #         # lnli2 <- -log((gamma(th+dni)/(factorial(dni)*gamma(th)))*(th/(th+sumnbcount))^th*(sumnbcount/(th+sumnbcount))^dni)
                                #         
                                #         # (gamma(r+k)/(factorial(k)*gamma(r)))*(r/(r+m))^r*(m/(r+m))^k})
                                #         # dnbinom(k, size=r, mu=m, log=FALSE)})
                                #         
                                #         
                                #         
                                #       } else {
                                #         nbjour <- data$ordinal2[i] - data$ordinal[i] + 1
                                #         # nbcount <- daily_count((1:nbjour) + data$ordinal[i] - 1, xparec, print = FALSE)
                                #         # Le premier jour de la série est 0 et donc il est à l'indice 1
                                #         nbcount <- nbperday[((1:nbjour) + data$ordinal[i] - 1)+1]
                                #         sumnbcount <- sum(nbcount)
                                #         lnli2 <- -dSnbinom(data$nombre[i], size = th, mu = nbcount, log = TRUE, tol=tol)
                                #       }
                                #       datatot[[k]]$LnL[i] <- lnli2
                                #       datatot[[k]]$Modeled[i] <- sumnbcount
                                #     }
                                #   } else {
                                #     for (i in 1:nrow(data)) {
                                #       if ((data$nombre[i] != 0) || zero_counts) {
                                #         if (is.na(data$Date2[i])) {
                                #           
                                #           # sumnbcount <- daily_count(data$ordinal[i], xparec, print = FALSE, zero = pt$zero)
                                #           # Le premier jour de la série est 0 et donc il est à l'indice 1
                                #           sumnbcount <- nbperday[data$ordinal[i]+1]
                                #           
                                #           # if (!is.null(pt$add.cofactors)) {
                                #           #   sumnbcount <- sumnbcount + sum(pt$cofactors[pt$cofactors$Date == 
                                #           #                                                 data$Date[i], pt$add.cofactor] * 
                                #           #                                    xparec[pt$add.cofactor])
                                #           #   if (sumnbcount <= pt$zero) 
                                #           #     sumnbcount <- pt$zero
                                #           # }
                                #           lnli2 <- -log(dnbinom(data$nombre[i], size = th, mu = sumnbcount, log = FALSE)/
                                #                            (1 - dnbinom(0, size = th, mu = sumnbcount,  log = FALSE)))
                                #           # dni <- data$nombre[i]
                                #           # lnli2 <- -log((gamma(th+dni)/(factorial(dni)*gamma(th)))*(th/(th+sumnbcount))^th*(sumnbcount/(th+sumnbcount))^dni*(1-((th/(th+sumnbcount))^th))^-1)
                                #           
                                #           # # After Arrabal CT, dos Santos Silva KP, Bandeira LN (2014) Zero-truncated negative binomial applied to nonlinear data. JP Journal of Biostatistics 11: 55-67 
                                #           # system.time({for (i in 1:100000) (gamma(r+k)/(factorial(k)*gamma(r)))*(r/(r+m))^r*(m/(r+m))^k*(1-((r/(r+m))^r))^-1})
                                #           # # After Girondot M (2010) Editorial: The zero counts. Marine Turtle Newsletter 129: 5-6 
                                #           # system.time({for (i in 1:100000) dnbinom(k, size=r, mu=m, log=FALSE)/(1-dnbinom(0, size=r, mu=m, log=FALSE))})
                                #           # 
                                #           
                                #         } else {
                                #           nbjour <- data$ordinal2[i] - data$ordinal[i] + 1
                                #           # nbcount <- daily_count((1:nbjour) + data$ordinal[i] - 1, xparec, print = FALSE)
                                #           # Le premier jour de la série est 0 et donc il est à l'indice 1
                                #           nbcount <- nbperday[((1:nbjour) + data$ordinal[i] - 1)+1]
                                #           sumnbcount <- sum(nbcount)
                                #           
                                #           # Est ce que les 0 ont été comptés
                                #           lnli2 <- -log(dSnbinom(data$nombre[i], size = th, mu = nbcount, log = FALSE, tol=tol)/
                                #                           (1 - dSnbinom(0, size = th, mu = nbcount, log = FALSE, tol=tol)))
                                #           
                                #           
                                #         }
                                #       } else {
                                #         lnli2 <- NA
                                #         sumnbcount <- NA
                                #       }
                                #       datatot[[k]]$LnL[i] <- lnli2
                                #       datatot[[k]]$Modeled[i] <- sumnbcount
                                #     }
                                #   }
                                # } else 
                                #   {
                                if (zero_counts) {
                                  
                                  for (i in 1:nrow(data)) {
                                    
                                    if (is.na(data$Date2[i])) {
                                      
                                      # sumnbcount <- daily_count(data$ordinal[i], xparec, print = FALSE, zero = pt$zero)
                                      # Le premier jour de la série est 0 et donc il est à l'indice 1
                                      sumnbcount <- nbperday[data$ordinal[i]+1]
                                      
                                      lnli2 <- -dnbinom(data$nombre[i], size = th, 
                                                        mu = sumnbcount, log = TRUE)
                                    } else {
                                      nbjour <- data$ordinal2[i] - data$ordinal[i] + 1
                                      # nbcount <- daily_count((1:nbjour) + data$ordinal[i] - 1, xparec, print = FALSE)
                                      # Le premier jour de la série est 0 et donc il est à l'indice 1
                                      nbcount <- nbperday[((1:nbjour) + data$ordinal[i] - 1)+1]
                                      sumnbcount <- sum(nbcount)
                                      lnli2 <- -dSnbinom(x=data$nombre[i], size = th, mu = nbcount, log = TRUE, tol=tol)
                                    }
                                    datatot[[k]]$LnL[i] <- lnli2
                                    datatot[[k]]$Modeled[i] <- sumnbcount
                                  }
                                } else {
                                  for (i in 1:nrow(data)) {
                                    if ((data$nombre[i] != 0) || zero_counts) {
                                      if (is.na(data$Date2[i])) {
                                        
                                        # sumnbcount <- daily_count(data$ordinal[i], xparec, print = FALSE, zero = pt$zero)
                                        # Le premier jour de la série est 0 et donc il est à l'indice 1
                                        sumnbcount <- nbperday[data$ordinal[i]+1]
                                        
                                        # if (!is.null(pt$add.cofactors)) {
                                        #   sumnbcount <- sumnbcount + sum(pt$cofactors[pt$cofactors$Date == 
                                        #                                                 data$Date[i], pt$add.cofactor] * 
                                        #                                    xparec[pt$add.cofactor])
                                        #   if (sumnbcount <= pt$zero) 
                                        #     sumnbcount <- pt$zero
                                        # }
                                        lnli2 <- -log(dnbinom(data$nombre[i], size = th, mu = sumnbcount, log = FALSE)/
                                                        (1 - dnbinom(0, size = th, mu = sumnbcount,  log = FALSE)))
                                      } else {
                                        nbjour <- data$ordinal2[i] - data$ordinal[i] + 1
                                        # nbcount <- daily_count((1:nbjour) + data$ordinal[i] - 1, xparec, print = FALSE)
                                        # Le premier jour de la série est 0 et donc il est à l'indice 1
                                        nbcount <- nbperday[((1:nbjour) + data$ordinal[i] - 1)+1]
                                        sumnbcount <- sum(nbcount)
                                        
                                        # Est ce que les 0 ont été comptés
                                        lnli2 <- -log(dSnbinom(data$nombre[i], size = th, mu = nbcount, log = FALSE, tol=tol)/
                                                        (1 - dSnbinom(0, size = th, mu = nbcount, log = FALSE, tol=tol)))
                                        
                                        
                                      }
                                    } else {
                                      lnli2 <- NA
                                      sumnbcount <- NA
                                    }
                                    datatot[[k]]$LnL[i] <- lnli2
                                    datatot[[k]]$Modeled[i] <- sumnbcount
                                  }
                                }
                                # }
                                datatot[[k]]$LnL[is.infinite(datatot[[k]]$LnL)] <- -log(1e-50)
                                return(datatot[[k]])
                                # }
                                
                                
                                
                              }, 
                              clusterExport = list(varlist=c("xpar", "datatot", "tol", 
                                                             "daily_count", "format_par", "pt"), 
                                                   envir=environment()), 
                              clusterEvalQ=list(expr=expression(library(phenology))
                              )
  ), silent = FALSE
  )
  
  if (any(class(rg) == "try-error")) {
    save(x, file = "x.Rdata")
    save(pt, file = "pt.Rdata")
    stop("Error during likelihood estimation; look at x.Rdata and pt.Rdata")
  }
  
  sum <- sum(sapply(X = rg, function(x) sum(x$LnL, na.rm = TRUE)))
  
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
