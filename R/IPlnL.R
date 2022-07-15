# IPlnL calculate the -log likelihood of data within a model.
# @title Calculate the -log likelihood of data within a model.
# @author Marc Girondot
# @return Return the -log likelihood of data within a model.\cr
# @param x A named vector of parameters.
# @param data Vector of observed interneting period.
# @param fixed.parameters Parameters that are fixed.
# @param verbose If TRUE, show the values of parameters at each step
# @param parallel If TRUE, parallel computing is used.
# @description Calculate the -log likelihood of data within a model.\cr
# @examples
# \dontrun{
# library(phenology)
# # Example
# data <- structure(c(`0` = 0, `1` = 47, `2` = 15, `3` = 6, `4` = 5, `5` = 4, 
#             `6` = 2, `7` = 5, `8` = 57, `9` = 203, `10` = 205, `11` = 103, 
#             `12` = 35, `13` = 24, `14` = 12, `15` = 10, `16` = 13, `17` = 49, 
#             `18` = 86, `19` = 107, `20` = 111, `21` = 73, `22` = 47, `23` = 30, 
#             `24` = 19, `25` = 17, `26` = 33, `27` = 48, `28` = 77, `29` = 83, 
#             `30` = 65, `31` = 37, `32` = 27, `33` = 23, `34` = 24, `35` = 22, 
#             `36` = 41, `37` = 42, `38` = 44, `39` = 33, `40` = 39, `41` = 24, 
#             `42` = 18, `43` = 18, `44` = 22, `45` = 22, `46` = 19, `47` = 24, 
#             `48` = 28, `49` = 17, `50` = 18, `51` = 19, `52` = 17, `53` = 4, 
#             `54` = 12, `55` = 9, `56` = 6, `57` = 11, `58` = 7, `59` = 11, 
#             `60` = 12, `61` = 5, `62` = 4, `63` = 6, `64` = 11, `65` = 5, 
#             `66` = 6, `67` = 7, `68` = 3, `69` = 2, `70` = 1, `71` = 3, `72` = 2, 
#             `73` = 1, `74` = 2, `75` = 0, `76` = 0, `77` = 3, `78` = 1, `79` = 0, 
#             `80` = 2, `81` = 0, `82` = 0, `83` = 1), 
#             Year = 1994, Species = "Dermochelys coriacea", 
#             location = "Yalimapo beach, French Guiana", class = "IP")
#             
# plot(data)
# 
# par <- c(ECF.2 = 0.042846438665386795, 
#          ECF.3 = 1.7321701716492339, 
#          ECF.4 = 2.2556217291284013, 
#          ECF.5 = 1.8117834932408947, 
#          ECF.6 = 1.3683652806822406, 
#          ECF.7 = 1.8844116917541645, 
#          ECF.8 = 0.72968863330728562, 
#          ECF.9 = 0.4691544793537184, 
#          ECF.10 = 0.058246650895644972, 
#          ECF.11 = 0.046330820923579881, 
#          ECF.12 = 0.0037895588708291416, 
#          ECF.13 = 0.0021687557575551016, 
#          ECF.14 = 0.00054829083259059011, 
#          meanIP = 9.9168825678697452, 
#          sdIP = 0.10142491419785879, 
#          minIP = 6.7979420131212036, 
#          pAbort = 2.4155971102805625, 
#          meanAbort = 1.9801272613187948, 
#          sdAbort = 0.5966756607139303, 
#          pCapture = -1.3270435936394733)
# 
# model <- IPModel(c(par, N=10000))
# 
# plot(model)
#      
# getFromNamespace(".IPlnL", ns="phenology")(par, fixed.parameters=c(N=10000), data=data)
# 
# }
# 
# @export

.IPlnL <- function (x, fixed.parameters, data, zero = 1e-10, verbose = FALSE, 
                    parallel = TRUE) 
{
  if (verbose) d(x)
  c <- IPModel(c(x, fixed.parameters), parallel = parallel)$cumuld
  if (is.matrix(data)) {
    if (ncol(data) > length(c)) {
      c <- c(c, rep(0, ncol(data) - length(c)))
    }        else {
      data <- cbind(data, matrix(0, ncol = length(c) - 
                                   ncol(data), nrow = nrow(data)))
    }
    pb <- c/sum(c)
    pb[pb == 0] <- zero
    if (any(is.na(pb))) {
      pb <- c(1, rep(zero, ncol(data)-1))
    }
    lnL <- 0
    for (r in 1:nrow(data)) {
      lnL <- lnL - dmultinom(data[r, ], prob = pb, log = TRUE)
    }
    return(lnL)
  }    else {
    if (any(is.na(c))) {
      pb <- c(1, rep(zero, length(data)-1))
      return(-dmultinom(data, prob = pb, log = TRUE))
    } else {
      if (length(data) > length(c)) {
        c <- c(c, rep(0, length(data) - length(c)))
      }        else {
        data <- c(data, rep(0, length(c) - length(data)))
      }
      pb <- c/sum(c)
      pb[pb == 0] <- zero
      
      return(-dmultinom(data, prob = pb, log = TRUE))
    }
  }
}
