.Tagloss_Lind <- function(individu, dfq, progressbar=FALSE) {
  
  p1 <- dfq$p1
  p2 <- dfq$p2
  pR1 <- dfq$pR1
  pR2 <- dfq$pR2
  pL1 <- dfq$pL1
  pL2 <- dfq$pL2
  Q1 <- dfq$Q1
  Q2 <- dfq$Q2
  Q3 <- dfq$Q3
  Q4 <- dfq$Q4
  Q5 <- dfq$Q5
  Q6 <- dfq$Q6
  Q7 <- dfq$Q7
  Q8 <- dfq$Q8
  LC_Q1 <- dfq$LC_Q1
  
  
  if (progressbar) {
    if (is.element('progress', installed.packages()[,1])) {
      # library("progress")
      pb <- getFromNamespace("progress_bar", ns="progress")$new(
        format = "  completion [:bar] :percent eta: :eta",
        total = nrow(individu), clear = FALSE)
      libp <- TRUE
    } else {
      libp <- FALSE
      pb <- txtProgressBar(min=0, max=nrow(individu), style=3)
    }
  }
  
  kl <- 0
  for (klj in 1:nrow(individu)) {
    
    if (progressbar) {
      if (libp) pb$tick() else setTxtProgressBar(pb, klj)
    }
    
    
    NLR_LR <- individu[klj, "NLR_LR"]
    NLR_L0 <- individu[klj, "NLR_L0"]
    NLR_0R <- individu[klj, "NLR_0R"]
    NL0_L0 <- individu[klj, "NL0_L0"]
    N0R_0R <- individu[klj, "N0R_0R"]
    NL0_00 <- individu[klj, "NL0_00"]
    N0R_00 <- individu[klj, "N0R_00"]
    NLR_00 <- individu[klj, "NLR_00"]
    
    N22 <- individu[klj, "N22"]
    N21 <- individu[klj, "N21"]
    N11 <- individu[klj, "N11"]
    N20 <- individu[klj, "N20"]
    N10 <- individu[klj, "N10"]  
    
    if ((!is.na(pR2[1]) & !is.na(pR1[1]) & !is.na(pL2[1]) & !is.na(pL1[1]))) {
      
      L_NLR_LR <- L_NLR_L0 <- L_NLR_0R <- L_NL0_L0 <- L_N0R_0R <- L_NL0_00 <- L_N0R_00 <- L_NLR_00 <- 1
      
      # NLR_LR
      if (!is.na(NLR_LR)) {
        L_NLR_LR <- prod(Q1[1:NLR_LR])
        } else {
          NLR_LR <- 0
        }
      # NLR_L0
      # k de NLR_LR + 1 à NLR_LR + NLR_L0
      if ((!is.na(NLR_L0) & !is.na(NLR_LR))) {
        L_NLR_L0_matrix <- matrix(data = c((NLR_LR+1):(NLR_LR+NLR_L0), rep(NA, 4*NLR_L0)), nrow=5, 
                                  dimnames=list(c("k", "Q1", "Q2", "Q5", "product")), byrow = TRUE)
        L_NLR_L0_matrix["Q2", ] <- Q2[L_NLR_L0_matrix["k", ]]
        L_NLR_L0_matrix["Q1", ] <- apply(X = L_NLR_L0_matrix, MARGIN=2, 
                                         FUN=function(x) prod(Q1[(NLR_LR + 1):(x["k"]-1)]))
        L_NLR_L0_matrix["Q5", ] <- apply(X = L_NLR_L0_matrix, MARGIN=2, 
                                         FUN=function(x) prod(Q5[(x["k"] + 1):(NLR_LR+NLR_L0)]))
        L_NLR_L0_matrix["Q5", ] <- ifelse(is.na(L_NLR_L0_matrix["Q5", ]), 1, L_NLR_L0_matrix["Q5", ])
        L_NLR_L0_matrix["product", ] <- apply(X= L_NLR_L0_matrix, MARGIN=2, 
                                              FUN = function(x) x["Q1"] * x["Q2"] * x["Q5"])
        L_NLR_L0 <- mean(L_NLR_L0_matrix["product", ])
      } else {
        NLR_L0 <- 0
      }
      # NLR_0R
      # k de NLR_LR + 1 à NLR_LR + NLR_0R
      if ((!is.na(NLR_0R) & !is.na(NLR_LR))) {
        L_NLR_0R_matrix <- matrix(data = c((NLR_LR+1):(NLR_LR+NLR_0R), rep(NA, 4*NLR_0R)), nrow=5, 
                                  dimnames=list(c("k", "Q1", "Q3", "Q7", "product")), byrow = TRUE)
        L_NLR_0R_matrix["Q3", ] <- Q3[L_NLR_0R_matrix["k", ]]
        L_NLR_0R_matrix["Q1", ] <- apply(X = L_NLR_0R_matrix, MARGIN=2, 
                                         FUN=function(x) prod(Q1[(NLR_LR + 1):(x["k"]-1)]))
        L_NLR_0R_matrix["Q7", ] <- apply(X = L_NLR_0R_matrix, MARGIN=2, 
                                         FUN=function(x) prod(Q7[(x["k"] + 1):(NLR_LR+NLR_0R)]))
        L_NLR_0R_matrix["Q7", ] <- ifelse(is.na(L_NLR_0R_matrix["Q7", ]), 1, L_NLR_0R_matrix["Q7", ])
        L_NLR_0R_matrix["product", ] <- apply(X= L_NLR_0R_matrix, MARGIN=2, 
                                              FUN = function(x) x["Q1"] * x["Q3"] * x["Q7"])
        L_NLR_0R <- mean(L_NLR_0R_matrix["product", ])
      } else {
        NLR_0R <- 0
      }
      # NL0_L0
      if ((!is.na(NLR_LR) & !is.na(NLR_L0) & !is.na(NL0_L0))) {
        L_NL0_L0 <- prod(Q5[(NLR_LR + NLR_L0 + 1):(NLR_LR + NLR_L0 + NL0_L0)])
      } else {
        NL0_L0 <- 0
      }
      # N0R_0R
      if ((!is.na(NLR_LR) & !is.na(NLR_0R) & !is.na(N0R_0R))) {
        L_N0R_0R <- prod(Q7[(NLR_LR + NLR_0R + 1):(NLR_LR + NLR_0R + N0R_0R)])
      } else {
        N0R_0R <- 0
      }
      # NL0_00
      if (!is.na(NLR_LR) & !is.na(NLR_L0) & !is.na(NL0_L0)  & !is.na(NL0_00)) {
        L_NL0_00_matrix <- matrix(data = c((NLR_LR + NLR_L0 + NL0_L0 + 1):(NLR_LR + NLR_L0 + NL0_L0 + NL0_00), rep(NA, 3*NL0_00)), 
                                  nrow=4, 
                                  dimnames=list(c("k", "Q5", "Q6", "product")), byrow = TRUE)
        L_NL0_00_matrix["Q6", ] <- Q6[L_NL0_00_matrix["k", ]]
        L_NL0_00_matrix["Q5", ] <- apply(X = L_NL0_00_matrix, MARGIN = 2, 
                                         FUN = function(x) prod(Q5[(NLR_LR + NLR_L0 + NL0_L0 + 1):(x["k"]-1)]))
        L_NL0_00_matrix["product", ] <- apply(X = L_NL0_00_matrix, MARGIN = 2, 
                                              FUN = function(x) x["Q6"] * x["Q5"])
        L_NL0_00 <- mean(L_NL0_00_matrix["product", ])
      } else {
        NL0_00 <- 0
      }
      # N0R_00
      if (!is.na(NLR_LR) & !is.na(NLR_0R) & !is.na(N0R_0R) & !is.na(N0R_00)) {
        L_N0R_00_matrix <- matrix(data = c((NLR_LR + NLR_0R + N0R_0R + 1):(NLR_LR + NLR_0R + N0R_0R + N0R_00), rep(NA, 3*N0R_00)), 
                                  nrow=4, 
                                  dimnames=list(c("k", "Q7", "Q8", "product")), byrow = TRUE)
        L_N0R_00_matrix["Q8", ] <- Q8[L_N0R_00_matrix["k", ]]
        L_N0R_00_matrix["Q7", ] <- apply(X = L_N0R_00_matrix, MARGIN = 2, 
                                         FUN = function(x) prod(Q7[(NLR_LR + NLR_0R + N0R_0R + 1):(x["k"]-1)]))
        L_N0R_00_matrix["product", ] <- apply(X = L_N0R_00_matrix, MARGIN = 2, 
                                              FUN = function(x) x["Q7"] * x["Q8"])
        L_N0R_00 <- mean(L_N0R_00_matrix["product", ])
      } else {
        N0R_00 <- 0
      }
      # NLR_00
      if (!is.na(NLR_00)) {
        L_NLR_00_matrix <- matrix(data = rep(0, NLR_00^2), ncol=NLR_00)
        
        for (JourPerteD in 1:NLR_00) {
          for (JourPerteG in 1:NLR_00) {
            
            if (JourPerteD == JourPerteG) {
              P <- ifelse(JourPerteG ==1, 1, prod(Q1[1:(JourPerteG-1)])) * Q4[JourPerteG]
            } else {
              if (JourPerteD < JourPerteG) {
                P <- ifelse(JourPerteD == 1, 1, prod(Q1[1:(JourPerteD-1)]))* Q2[JourPerteD] * prod(Q5[(JourPerteD+1):(JourPerteG-1)]) * Q6[JourPerteG]
              } else {
                P <- ifelse(JourPerteG == 1, 1, prod(Q1[1:(JourPerteG-1)]))* Q3[JourPerteG] * prod(Q7[(JourPerteG+1):(JourPerteD-1)]) * Q8[JourPerteG]
              }
            }
            
            L_NLR_00_matrix[JourPerteG, JourPerteD] <- P
          }
        }
        
        L_NLR_00 <- mean(L_NLR_00_matrix)
      }
      
      L <- log(L_NLR_LR) + log(L_NLR_L0) + log(L_NLR_0R) + log(L_NL0_00) + log(L_N0R_00) + log(L_NLR_00)
      
    } else {
      if ((!is.na(p2[1]) & !is.na(p1[1]))) {
        
        # N22 signifie que l'individu a été vu avec 2 bagues du jour 1 au jour N22 inclus
        # Si N22=0 c'est que le lendemain du jour bagué, il a été revu avec 1 bague
        if ((!is.na(N22) & (N22 != 0))) {
          # Q22 <- prod(Q1_[1:N22]) 
          LQ22 <- LC_Q1[N22]
        } else {
          LQ22 <- 0
        }
        
        # N21 C'est le nombre de jours entre le vu avec 2 bagues et avec une seule
        # Si N21=0
        if (!is.na(N22) & !is.na(N21) & (N21 != 0)) {
          Q21_matrix <- matrix(data = c((N22+1):(N22+N21), rep(NA, N21*4)), nrow=5, 
                               dimnames=list(c("k", "Q1", "Q3", "Q4", "product"), NULL), 
                               byrow = TRUE)
          Q21_matrix["Q1", ] <- apply(X = Q21_matrix, MARGIN = 2, 
                                      FUN = function(x) prod(Q1[(N22+1):(x["k"]-1)]))
          Q21_matrix["Q1", 1] <- 1
          Q21_matrix["Q3", ] <- Q3[Q21_matrix["k", ]]
          Q21_matrix["Q4", ] <- apply(X = Q21_matrix, MARGIN = 2, 
                                      FUN = function(x) prod(Q4[(x["k"]+1):(N21 + N22)]))
          Q21_matrix["Q4", N21] <- 1
          #      Q21_matrix["Q4", ] <- ifelse(is.na(Q21_matrix["Q4", ]), 1, Q21_matrix["Q4", ])
          Q21_matrix["product", ] <- Q21_matrix["Q1", ] * Q21_matrix["Q3", ] * Q21_matrix["Q4", ]
          Q21 <- mean(Q21_matrix["product", ])
        } else {
          Q21 <- 1
        }
        
        if (!is.na(N22) & !is.na(N21) & !is.na(N11) & (N11 != 0)) {
          Q11 <- prod(Q4[(N22 + N21 + 1):(N22 + N21 + N11)])
        } else {
          Q11 <- 1
        }
        
        if (!is.na(N22) & !is.na(N21) & !is.na(N11) & !is.na(N10) & (N10 != 0)) {
          Q10_matrix <- matrix(data = c((N22+N21+N11+1):(N22+N21+N11+N10), rep(NA, N10*3)), nrow=4, 
                               dimnames=list(c("k", "Q4", "Q5", "product"), NULL), 
                               byrow = TRUE)
          Q10_matrix["Q4", ] <- apply(X = Q10_matrix, MARGIN = 2, 
                                      FUN = function(x) prod(Q4[(N22+N21+N11+1):(x["k"]-1)]))
          Q10_matrix["Q4", 1] <- 1
          Q10_matrix["Q5", ] <- Q5[Q10_matrix["k", ]]
          Q10_matrix["product", ] <- Q10_matrix["Q4", ] * Q10_matrix["Q5", ]
          Q10 <- mean(Q10_matrix["product", ])
        } else {
          Q10 <- 1
        }
        
        if (!is.na(N20) & (N20 != 0)) {
          Q20_matrix <- matrix(data = rep(0, N20^2), ncol=N20)
          
          for (JourPerteD in 1:N20) {
            for (JourPerteG in 1:N20) {
              
              if (JourPerteD == JourPerteG) {
                P <- ifelse(JourPerteG ==1, 1, prod(Q1[1:(JourPerteG-1)])) * Q2[JourPerteG]
              } else {
                if (JourPerteD < JourPerteG) {
                  P <- ifelse(JourPerteD == 1, 1, prod(Q1[1:(JourPerteD-1)]))* Q3[JourPerteD] * ifelse((JourPerteG-1)-(JourPerteD+1)==0, 1, prod(Q4[(JourPerteD+1):(JourPerteG-1)])) * Q5[JourPerteG]
                } else {
                  P <- ifelse(JourPerteG == 1, 1, prod(Q1[1:(JourPerteG-1)]))* Q3[JourPerteG] * ifelse((JourPerteD-1)-(JourPerteG+1)==0, 1, prod(Q4[(JourPerteG+1):(JourPerteD-1)])) * Q5[JourPerteG]
                }
              }
              
              Q20_matrix[JourPerteG, JourPerteD] <- P
            }
          }
          
          Q20 <- mean(Q20_matrix)
          
        } else {
          Q20 <- 1
        }
        
        L <- LQ22 + log(Q21) + log(Q11) + log(Q10) + log(Q20)
      }
  }
    # print(paste(kl, "   ", L))
    kl <- kl + L 
  }
  return(kl)
}