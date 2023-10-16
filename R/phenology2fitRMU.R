#' phenology2fitRMU is used to prepare output of phenology to be used by fitRMU()
#' @title Create the data to be used with fitRMU() for a summary(phenology).
#' @author Marc Girondot
#' @return Return a list with elements ready to be used with fitRMU().
#' @param phenologyout A result of plot.phenology() or summary()
#' @param col.mean Name of the column to be used as mean value. Can be a vector.
#' @param col.var Name of the column to be used as variance value. Can be a vector.
#' @param rookeries.names.grep The pattern to return the rookery name from names of timeseries.
#' @param years.grep The pattern to return the year from names of timeseries.
#' @param limit.cv Remove data with higher coefficient of variation than the limit.
#' @param limit.sd Remove data with higher standard deviation than the limit.
#' @param density What density should be used. Can be dnorm or dgamma. Can be a vector.
#' @description This function takes the result of plot.phenology() and generates the information 
#' to be used with fitRMU().\cr
#' The value of density can be dnorm or dgamma. dnorm is better if ML results are used and 
#' dgamma is for MCMC.\cr
#' Here are some example of regular expressions (regex) in grep:\cr
#' If format of timeseries is beachname-2005: \cr
#' rookeries.names.grep="(.+)-.+"\cr
#' years.grep=".+-(.+)$"\cr
#' Examples:\cr
#' gsub("(.+)-.+", "\\\\1", "beachname-2005"); gsub(".+-(.+)$", "\\\\1", "beachname-2005")\cr
#' If format of timeseries is beachname-2005-2006: \cr
#' rookeries.names.grep="(.+)-.+-.+"\cr
#' years.grep=".+-([0-9][0-9][0-9][0-9])-.+$"\cr
#' Examples:\cr
#' gsub("(.+)-.+-.+", "\\\\1", "beachname-2005-2006"); 
#' gsub(".+-([0-9][0-9][0-9][0-9])-.+$", "\\\\1", "beachname-2005-2006")\cr
#' If format of timeseries is beachname-20052006:\cr
#' rookeries.names.grep="(.+)-.+"\cr
#' years.grep=".+-([0-9][0-9][0-9][0-9])([0-9][0-9][0-9][0-9])$"\cr
#' Examples:\cr
#' gsub("(.+)-.+", "\\\\1", "beachname-20052006"); 
#' gsub(".+-([0-9][0-9][0-9][0-9])([0-9][0-9][0-9][0-9])$", "\\\\1", "beachname-20052006")\cr
#' The return is a list with these elements:\cr
#' RMU.data, years.byrow, colname.year, and RMU.names.\cr
#' If density is a vector, the density used is linked to the rank of the timeseries in phenologyout.
#' @family Fill gaps for RMU
#' @family Phenology model
#' @examples
#' \dontrun{
#' library("phenology")
#' }
#' @export

phenology2fitRMU <- function(phenologyout=stop("A result obtained from summary(phenology)"), 
                             col.mean="with_obs_Mean_ML", 
                             col.var="with_obs_Var_ML", 
                             rookeries.names.grep="(.+)-.+",
                             years.grep=".+-(.+)$", 
                             limit.cv=+Inf,
                             limit.sd=+Inf,
                             density="dgamma") {
  synthesis <- phenologyout$synthesis
  density <- rep(density, nrow(synthesis))[1:nrow(synthesis)]
  
  if (!is.null(col.var)) {
    synthesis <- cbind(synthesis, sd=sqrt(synthesis[, col.var]))
  } else {
    synthesis <- cbind(synthesis, sd=rep(0, nrow(synthesis)))
  }
  synthesis <- cbind(synthesis, cv=synthesis[, "sd"]/synthesis[, col.mean])
  synthesis <- synthesis[synthesis$cv<limit.cv, ]
  synthesis <- synthesis[synthesis$sd<limit.sd, ]
  col.sd <- "sd"
  col.density <- "density"
  
  # 17/5/2020: je laisse plus de solutions à l'utilisateur 
  
  liste_motus <- gsub(rookeries.names.grep, "\\1", synthesis$series)
  liste_years <- gsub(years.grep, "\\1", synthesis$series)
  
  motus  <- unique(liste_motus)
  years <- unique(liste_years)
  years <- years[order(as.numeric(years))]
  
  df.out <- data.frame(Year=as.numeric(years), stringsAsFactors = FALSE, row.names = years)
  
  for (i in seq_along(motus)) {
    dfi <- data.frame(motu.mean=rep(NA, nrow(df.out)), 
                      motu.SD=rep(NA, nrow(df.out)), 
                      motu.density=rep(NA, nrow(df.out)), stringsAsFactors = FALSE
    )
    colnames(dfi) <- c(paste0(motus[i], ".mean"), 
                       paste0(motus[i], ".sd"), 
                       paste0(motus[i], ".density"))
    df.out <- cbind(df.out, dfi, stringsAsFactors = FALSE)
  }
  # Maintenant je remplis
  
  for (n in 1:nrow(synthesis)) {
    df.out[liste_years[n], paste0(liste_motus[n], ".mean")] <- synthesis[n, col.mean]
    df.out[liste_years[n], paste0(liste_motus[n], ".sd")] <- synthesis[n, col.sd]
    df.out[liste_years[n], paste0(liste_motus[n], ".density")] <- density[n]
  }
  
  dfRMU.names <- data.frame(mean=colnames(df.out)[grepl(".mean$", colnames(df.out))], 
                            se=colnames(df.out)[grepl(".sd$", colnames(df.out))], 
                            density=colnames(df.out)[grepl(".density$", colnames(df.out))], stringsAsFactors = FALSE)
  
  listout <- list(RMU.data=df.out, 
                  years.byrow=TRUE, 
                  colname.year="Year", 
                  RMU.names=dfRMU.names)
  
  return(listout)
}

