#' TableECFOCF formats a CMR dataset into a file that fitCF can use.
#' @title Format a CMR dataset into a file that fitCF can use.
#' @author Marc Girondot
#' @return Return a matrix with counts for all OCF and ECF combinations.
#' @param data CMR file.
#' @param columnID Name of the column in data for unique identifier of females.
#' @param columnDate Name of the column in data for morning date when female has been seen on the beach.
#' @param MinimumDaysBetween2Nest Number of minimum days between two nests.
#' @param MeanDaysBetween2Nest Average number of days between two nests.
#' @param MaxNests Maximum number of nests by a female.
#' @param date0 Date for the ordinal day 0.
#' @param length_season The total length of season based on groups of interclutch intervals.
#' @description This function formats a CMR dataset to a file that fitCF can use.\cr
#' If date0 is not null, a 3D TableECFOCF is generated.\cr
#' 3D table (ECF, OCF, period) has two attributes:\cr
#' - table with 4 elements:\cr
#' begin, end are the first and last elements with counts\cr
#' min and max are the first and last period where a nest could have been laid based on MaxNests value\cr
#' - characteristics with 5 elements:\cr
#' MinimumDaysBetween2Nest, MeanDaysBetween2Nest MaxNests, date0, length_season\cr
#' p parameter can be setup to +Inf until begin and after end 
#' @family Model of Clutch Frequency
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' data(MarineTurtles_2002)
#' ECFOCF_2002 <- TableECFOCF(MarineTurtles_2002)
#' plot(ECFOCF_2002)
#' ECFOCF_2002 <- TableECFOCF(MarineTurtles_2002, date0=as.Date("2002-01-01"))
#' plot(ECFOCF_2002, period=11)
#' }
#' @export

# Table ECF OCF ####

TableECFOCF <- function(data=stop("A dataframe with a column 'ID' and a column 'Date'"), 
                        columnID="ID", 
                        columnDate="Date", 
                        MinimumDaysBetween2Nest=7, 
                        MeanDaysBetween2Nest=9.8, 
                        MaxNests=15, 
                        date0=NULL, 
                        length_season=floor(365/MeanDaysBetween2Nest)+1) {
  
  nf <- levels(as.factor(data[, columnID]))
  
  ECFOCF_Obs <- data.frame(ID=character(), 
                           OCF=numeric(), 
                           ECF=numeric(), 
                           period=numeric(), 
                           stringsAsFactors = FALSE)
  for (i in 1:length(nf)) {

    sb <- subset(x = data, subset=data[, columnID]==nf[i], 
                 select = columnDate, drop=TRUE)
    if (!is.null(date0)) {
      period <- floor(as.numeric(sb[order(sb)][1] - date0 + 1) / MeanDaysBetween2Nest) + 1
    } else {
      period <- 1
    }
    
    if (length(sb) == 1) {
      OCF=1
      ECF=1
    } else {
      delta <- as.numeric(diff(sb[order(sb)]))
      OCF = length(delta[delta>=MinimumDaysBetween2Nest])+1
      ECF = floor((sum(delta)/MeanDaysBetween2Nest)+1)
    }
    if (OCF>ECF) ECF <- OCF
    ECFOCF_Obs <- rbind(ECFOCF_Obs, 
                        data.frame(ID=nf[i], 
                                   OCF=OCF,
                                   ECF=ECF, 
                                   period=period, 
                                   stringsAsFactors = FALSE))
  }
  
  end <- max(ECFOCF_Obs$period)
  begin <- min(ECFOCF_Obs$period)
  
  
  if (!is.null(date0)) {
  OCFECF <- array(data = 0, dim=c(MaxNests+1, MaxNests+1, length_season+MaxNests+1), 
                  dimnames = list(paste0("OCF", 0:(MaxNests)), paste0("ECF", 0:(MaxNests)), 
                                  paste0("time", 1:(length_season+MaxNests+1))
                  )
  )
  } else {
    OCFECF <- array(data = 0, dim=c(MaxNests+1, MaxNests+1, 1), 
                    dimnames = list(paste0("OCF", 0:(MaxNests)), paste0("ECF", 0:(MaxNests)), 
                                    paste0("time", 1)
                    )
    )
  }
  
  for (i in 1:nrow(ECFOCF_Obs)) {
    OCFECF[ECFOCF_Obs[i, "OCF"]+1, ECFOCF_Obs[i, "ECF"] +1, ECFOCF_Obs[i, "period"]] <- 1 + OCFECF[ECFOCF_Obs[i, "OCF"]+1, ECFOCF_Obs[i, "ECF"] +1, ECFOCF_Obs[i, "period"]]
  }
  
 class(OCFECF) <- "TableECFOCF"
 
 if (end==1) {
   attributes(OCFECF) <- c(attributes(OCFECF), table=list(c(begin=1, end=1, min=1, max=1)))
 } else {
 # part of the table for which data are present
 # range <- which(sapply(1:(dim(OCFECF)[3]), function(dim3) {sum(OCFECF[, , dim3])!=0}))
 attributes(OCFECF) <- c(attributes(OCFECF), table=list(c(begin=begin, end=end, min=begin-MaxNests, max=end+MaxNests)))
 }
 
 attributes(OCFECF) <- c(attributes(OCFECF), characteristics=list(c(MinimumDaysBetween2Nest=MinimumDaysBetween2Nest, 
                                                                    MeanDaysBetween2Nest=MeanDaysBetween2Nest,
                                                                    MaxNests=MaxNests, 
                                                                    date0=date0, 
                                                                    length_season=length_season)))
  return(OCFECF)
}

