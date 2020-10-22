.LikelihoodRMU <- function (x, fixed.parameters, 
                            model.trend, 
                            model.SD, 
                            RMU.data, 
                            colname.year = NULL, 
                            RMU.names = NULL, 
                            index = NULL) {
  
  # index <- list(year = index.year, mean = index.mean, se = index.se, 
  #               density = index.density, 
  #               colnames = nabeach, nyear = nyear, nbeach = nbeach, maxL = maxL)
  
  
  x <- c(x, fixed.parameters)
  
  # print(x)
  
  if (is.null(index)) {
    nm <- colnames(RMU.data)
    index.year <- which(nm == colname.year)
    index.mean <- match(RMU.names$mean, nm)
    if (!is.null(RMU.names$se)) {
      index.se <- match(RMU.names$se, nm)
    } else {
      index.se <- NULL
    }
    if (!is.null(RMU.names$density)) {
      index.density <- match(RMU.names$density, nm)
    } else {
      index.density <- NULL
    }
    d <- RMU.data[, index.mean]
    nabeach <- colnames(d)
    nbeach <- length(nabeach)
    nyear <- dim(RMU.data)[1]
    nayear <- RMU.data[, index.year]
    maxL <- 1e+09
    index <- list(year = index.year, mean = index.mean, se = index.se, 
                  density=NULL, 
                  colnames = nabeach, nyear = nyear, nbeach = nbeach)
  } else {
    nabeach <- index$colnames
    nbeach <- length(nabeach)
    nyear <- dim(RMU.data)[1]
    maxL <- index$maxL
    index.density <- index$density
  }
  
  SD <- abs(x[(substr(names(x), 1, 3) == "SD_")])
  if ((length(SD) == 0) | (model.SD == "zero")) {
    SD <- rep(0, nbeach)
    names(SD) <-  paste0("SD_", nabeach)
  }    else {
    if (length(SD) == 1) {
      SD <- rep(SD, nbeach)
      names(SD) <-  paste0("SD_", nabeach)
   }
  }
  SD <- SD[paste0("SD_", nabeach)]
  
  aSD <- abs(x[(substr(names(x), 1, 4) == "aSD_")])
  if ((length(aSD) == 0) | (model.SD == "zero")) {
    aSD <- rep(0, nbeach)
    names(aSD) <-  paste0("aSD_", nabeach)
  }    else {
    if (length(aSD) == 1) {
      aSD <- rep(aSD, nbeach)
      names(aSD) <- paste0("aSD_", nabeach)
    }
  }
  aSD <- aSD[paste0("aSD_", nabeach)]
  
  dtaL_obs <- RMU.data[, index$mean]
  # 20/9/2019
  if (!is.null(index.density)) {
    dtaL_density <- RMU.data[, index.density]
  } else {
    dtaL_density <- NULL
  }
  
  if (!is.null(index$se)) {
    # Il n'y a pas de SE observé
    
    dtaL_SD <- sqrt(RMU.data[, index$se]^2 + (dtaL_obs * 
                                                matrix(rep(aSD, nyear), nrow = nyear, byrow = TRUE) + 
                                                matrix(rep(SD, nyear), nrow = nyear, byrow = TRUE))^2)
  }    else {
    dtaL_SD <- (dtaL_obs * 
                  matrix(rep(aSD, nyear), nrow = nyear, byrow = TRUE) + 
                  matrix(rep(SD, nyear), nrow = nyear, byrow = TRUE))
  }
  
  La0 <- x[paste0("a0_", nabeach[paste0("a0_", nabeach) %in% names(x)])]
  La1 <- x[paste0("a1_", nabeach[paste0("a1_", nabeach) %in% names(x)])]
  La2 <- x[paste0("a2_", nabeach[paste0("a2_", nabeach) %in% names(x)])]
  if (any(is.na(La2))) {
    La2 <- La0
    La2[] <- 0
    names(La2) <- gsub("^a0_", "a2_", names(La0))
  }
  if (any(is.na(La1))) {
    La1 <- La0
    La1[] <- 0
    names(La1) <- gsub("^a0_", "a1_", names(La0))
  }
  map <- matrix(rep(NA, nyear * (nbeach)), ncol = nbeach)
  
  # Cela va de 1 à x
  nye <- RMU.data[, colname.year]-min(RMU.data[, colname.year])+1
  rownames(map) <- names(nye)
  
  for (i in 1:nbeach) {
    map[, i] <- abs(La2[i]) * nye^2 + abs(La1[i]) * nye + abs(La0[i])
  }
  mapp <- matrix(rep(NA, nyear * nbeach), ncol = nbeach, dimnames = list(names(nye), 
                                                                         beach = nabeach))
  for (j in 1:nyear) {
    mapp[j, ] <- map[j, ] / sum(map[j, ])
  }
  if (model.trend == "year-specific") {
    Tot <- abs(x[paste0("T_", RMU.data[, index$year])])
    dtaL_theo <- matrix(rep(Tot, nbeach), ncol = nbeach, 
                        byrow = FALSE, dimnames = list(NULL, beach = nabeach))
  }
  if (model.trend == "constant") {
    Tot <- abs(x["T_"])
    dtaL_theo <- matrix(rep(Tot, nbeach * nyear), ncol = nbeach, 
                        byrow = FALSE, dimnames = list(NULL, beach = nabeach))
  }
  if (model.trend == "exponential") {
    Tot <- abs(x["T_"]) * exp(x["r"] * (nye))
    dtaL_theo <- matrix(rep(Tot, nbeach), ncol = nbeach, 
                        byrow = FALSE, dimnames = list(names(nye), beach = nabeach))
  }
  dtaL_obs <- dtaL_obs[, nabeach]
  dtaL_theo <- dtaL_theo * mapp
  valide <- !is.na(dtaL_obs)
  for (i in 1:nrow(dtaL_SD)) dtaL_SD[i, (dtaL_SD[i, ] == 0) & 
                                       (!is.na(dtaL_SD[i, ]))] <- 1e-10
  if (any(is.na(dtaL_theo))) {
    L <- maxL
  }    else {
    # C'est là où je dois faire du dnorm ou du dgamma
    
    if (any(dtaL_density == "dnorm", na.rm = TRUE) | (is.null(dtaL_density))) {
      # save(valide, dtaL_obs, dtaL_density, dtaL_SD, file="Atester.Rdata")
      
      if (is.null(dtaL_density)) {
        L.dnorm <- sum(-dnorm(x = dtaL_obs[valide], 
                              mean = dtaL_theo[valide], 
                              sd = dtaL_SD[valide], log = TRUE))
      } else {
        L.dnorm <- sum(-dnorm(x = dtaL_obs[valide & (dtaL_density == "dnorm")], 
                              mean = dtaL_theo[valide & (dtaL_density == "dnorm")], 
                              sd = dtaL_SD[valide & (dtaL_density == "dnorm")], log = TRUE))
      }
    } else {
      L.dnorm <- 0
    }
    
    if (any(dtaL_density[valide] == "dgamma", na.rm = TRUE)) {
      scale <- (dtaL_SD[valide & (dtaL_density == "dgamma")])^2 / dtaL_theo[valide & (dtaL_density == "dgamma")]
      shape <- (dtaL_theo[valide & (dtaL_density == "dgamma")] ^2) / ((dtaL_SD[valide & (dtaL_density == "dgamma")])^2)
      # rate <- 1/scale
      
      L.gamma <- sum(-dgamma(x = dtaL_obs[valide & (dtaL_density == "dgamma")], 
                             shape = shape, 
                             scale = scale, log = TRUE))
    } else {
      L.gamma <- 0
    }
    L <- L.dnorm + L.gamma
  }
  # print(L)
  return(L)
}

# je donne les N probabilités et j'ai les N-1 probabilités conditionelles 
# Ne sert plus
.inv.p.multinomial <- function(x) {
  p <- x[1]
  for (i in 2:(length(x)-1)) 
    p <- c(p, x[i]/(prod(1-p[1:(i-1)])))
  return(p)
}

# je donne les N-1 probabilités conditionelles et j'ai les N probabilités
.p.multinomial <- function(p) c(p, 1)*c(1, sapply(seq_along(p), function(i) prod(1-p[1:i])))

